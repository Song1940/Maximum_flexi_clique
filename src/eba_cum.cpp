#include "eba.h"
#include "fpa.h"
#include "mgpp.h"
#include "seed_heuristic.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <queue>
#include <utility>
#include <vector>
#include <cstdio>
#ifdef RANDOM_ORDER
#include <random>
#include <cstdlib>
#endif

// ============================================================
// eba_cum.cpp — cumulative connectivity-preserving branching EBA
// (exact maximum flexi-clique). A search node partitions vertices into
// S / Cr / Cun / D; CP-branching enumerates connected subgraphs with a
// saturation prune and an optional width-based pivot.
//
// CP-branching (§2.2): order the reachable candidates Cr = <v_0..v_{k-1}>.
// A node spawns k+1 children: for each exclusion child B_i the partial set
// S accumulates the prefix {v_0..v_{i-1}} while v_i alone is excluded; the
// include-all child B_k admits the whole ordering.  Children partition all
// targets by *first excluded candidate*, so enumeration is complete and
// duplicate-free.  Realised here as a single sweep that pushes one node into
// S per step, so the incremental degree machinery applies:
//
//   for idx, v_i in order:            // S already == nodeS ∪ {v_0..v_{idx-1}}
//       exclude v_i  -> recurse       // child B_i
//       restore v_i; push v_i -> S    // accumulate prefix for next sibling
//   recurse                           // include-all child B_k
//   pop the accumulated prefix
//
// Pruning — the paper's two rules (§5.2), applied per node in refine(), plus a
// size-bound termination check:
//   Rule 1 (degree):   drop candidates / prune the node by adjusted degree;
//   Rule 2 (diameter): drop candidates too far from / disconnected from S;
//   size bound:        prune once |S| reaches its size bound while infeasible.
// A global form of Rule 1 (degree-core reduction) is also run at preprocessing
// and is tightened whenever F' improves. Connectivity is automatic (every Cr
// node is adjacent to S), so it is never checked.
//
// Width-based pivot (§5.3, default ON; NO_PIVOT to ablate): a width-control candidate
// ordering, NOT a separate rule. It picks the least-slack node p and front-loads
// p's non-neighbours so the cumulative prefix drives p's in-S non-neighbour count
// dbar_S(p) up to its tolerance dbar_max(p); the moment p saturates the sibling
// loop is cut, bounding branching width. The saturation test it relies on,
//   dbar_max(v) = s_max(v) - 1 - floor(s_max(v)^τ),  s_max(v)=max(1,min(⌊(d_a+1)^{1/τ}⌋,|R|))
//   dbar_S(v)   = (|S|-1) - deg_in_S(v),
// (prune when s_max(v) < |S| or dbar_S(v) > dbar_max(v)) is sound on its own —
// S ⊆ H implies dbar_H(v) ≥ dbar_S(v) — and is checked in refine(). The pivot is
// a pure reordering and never changes the returned maximum size.
// ============================================================

#ifdef TRACE_HIST
long long g_hist[256];
#endif

#ifdef TRACE_WIDTH
// Branching-width histogram (Appendix, Thm. "conditional bound below 2^n"):
// g_width_hist[w] = number of CP-branching search nodes that explored exactly w
// children (w = explored exclusion children + the include-all child B_k, if any).
// Only candidate-branching nodes are counted (not leaves / pruned / the S=0 seed
// loop), since the width-w recurrence T(mu) <= sum_{j=1}^{w} T(mu-j) of Thm. C.x
// is stated for those nodes. Default OFF; does not alter the search itself.
long long g_width_hist[256];
#endif

namespace {

struct CumState {
    const Graph& G;
    double tau;

    std::vector<int> best;
    int best_size = 0;

    // Partition / scope state (indexed by vertex).
    std::vector<bool> in_S;       // v ∈ S
    std::vector<bool> in_D;       // v ∈ D (excluded)
    std::vector<bool> in_scope;   // v ∈ R = S ∪ Cr ∪ Cun
    std::vector<int>  local_adj_deg;  // d_a(v) = |adj[v] ∩ R|
    std::vector<int>  deg_in_S;       // |adj[v] ∩ S|

    // Rule 2 (diameter pruning): dist_S[u] = max_{v∈S} sp_{G[H]}(u,v),
    // the longest shortest-path from u to any S vertex within the induced
    // subgraph G[H]. 0 when S=∅; INT_MAX when u is disconnected from S in G[H].
    // Maintained incrementally (one BFS per pushToS, max-combined, undone on pop).
    // NOTE: distances are NOT re-raised when a later scopeRemove shrinks G[H], so a
    // stored finite dist_S[u] may underestimate the true sp in the reduced G[H].
    // This only ever weakens Rule 2 (an underestimate prunes less), never
    // over-prunes, so soundness holds; some disconnection prunes are just missed.
    std::vector<int>  dist_S;
    // BFS scratch reused via touched-list resets (bfs_d kept all-INT_MAX between calls).
    std::vector<int>  bfs_d, bfs_touched, bfs_queue;

    long long branches = 0;

#ifdef RANDOM_ORDER
    // RNG for the random-ordering pivot ablation. Reproducible per seed: the seed
    // defaults to 0xC0FFEE but can be overridden via the FLEXI_SEED env var, so the
    // ablation can measure variance across multiple random orderings.
    std::mt19937 rng{ []() -> unsigned {
        const char* s = std::getenv("FLEXI_SEED");
        return s ? static_cast<unsigned>(std::strtoul(s, nullptr, 10)) : 0xC0FFEEu;
    }() };
#endif

    CumState(const Graph& G, double tau)
        : G(G), tau(tau),
          in_S(G.n, false), in_D(G.n, false), in_scope(G.n, false),
          local_adj_deg(G.n, 0), deg_in_S(G.n, 0),
          dist_S(G.n, 0), bfs_d(G.n, INT_MAX) {}

    int theta() const { return computeTheta(best_size, tau); }

    // floor(x^{1/τ}), clamped to INT_MAX.
    // With small τ the exponent 1/τ is large (e.g. τ=0.2, x=128 → 128^5 ≈ 3.4e10),
    // which overflows int; casting an out-of-range double to int is UB. Every
    // caller either caps the result by scope_size (≤ n ≤ INT_MAX) or compares it
    // against |S| (≤ n), so saturating at INT_MAX is value-preserving here while
    // removing the UB.
    int floorRoot(int x) const {
        if (x <= 0) return 0;
        double r = std::floor(std::pow((double)x, 1.0 / tau));
        if (r >= (double)INT_MAX) return INT_MAX;
        return (int)r;
    }

    // S is always connected (CP invariant), so only the degree test is needed.
    bool feasibleS(const std::vector<int>& S) const {
        int k = (int)S.size();
        if (k == 0) return false;
        int thr = floorPow(k, tau);
        for (int v : S)
            if (deg_in_S[v] < thr) return false;
        return true;
    }

    // Remove v from scope R; record undo (decremented neighbours + removed v).
    void scopeRemove(int v, std::vector<std::pair<int,int>>& deg_undo,
                     std::vector<int>& removed) {
        in_scope[v] = false;
        removed.push_back(v);
        for (int u : G.adj[v])
            if (in_scope[u]) { deg_undo.push_back({u, local_adj_deg[u]}); --local_adj_deg[u]; }
    }

    // Undo scopeRemove/cascade: restore degrees (reverse) then scope membership.
    void undoScope(std::vector<int>& removed,
                   std::vector<std::pair<int,int>>& deg_undo,
                   bool reset_in_D = false) {
        for (auto it = deg_undo.rbegin(); it != deg_undo.rend(); ++it)
            local_adj_deg[it->first] = it->second;
        for (int v : removed) { in_scope[v] = true; if (reset_in_D) in_D[v] = false; }
    }

    // Move v into S (still in scope); update incremental deg_in_S.
    // deg_in_S[u] = |adj[u] ∩ S| is maintained for EVERY vertex u (not just S
    // members), because buildChildSets uses it to tell Cr (adjacent to S) from
    // Cun. v's own count is already accumulated by earlier pushes of its
    // neighbours, so only its neighbours are bumped here.
    void pushToS(int v, std::vector<int>& S) {
        in_S[v] = true;
        S.push_back(v);
        for (int u : G.adj[v]) ++deg_in_S[u];
    }
    void popFromS(int v, std::vector<int>& S) {
        in_S[v] = false;
        S.pop_back();
        for (int u : G.adj[v]) --deg_in_S[u];
    }

    // Non-neighbours of v inside S: |S| (minus 1 if v∈S) − |adj[v] ∩ S|.
    int nonDegInS(int v, int s_size) const {
        return (s_size - (in_S[v] ? 1 : 0)) - deg_in_S[v];
    }

    int dbarMax(int v, int scope_size) const {
        int da = local_adj_deg[v];
        int smax = std::max(1, std::min(floorRoot(da + 1), scope_size));
        return smax - 1 - floorPow(smax, tau);
    }
    int smaxOf(int v, int scope_size) const {
        int da = local_adj_deg[v];
        return std::max(1, std::min(floorRoot(da + 1), scope_size));
    }

    // Lemma (degree-diameter): a connected graph with min degree k and diameter L
    // has at least n(k,L) nodes. Long-long version of graph.h::nFunc to avoid the
    // (L/3)*(k-2) term overflowing int for large H.
    static long long nFuncLL(int k, int L) {
        if (L <= 0) return 1;
        if (k <= 0) return (long long)L + 1;
        if (L <= 2 || k == 1) return (long long)k + L;
        return (long long)k + L + 1 + (long long)(L / 3) * (k - 2);
    }

    // BFS from w within G[H] (in_scope vertices); update dist_S[u] = max(dist_S[u],
    // sp(w,u)) for every u in `universe` (INT_MAX when w cannot reach u in G[H]).
    // Records each raised value in dist_undo for exact restoration on pop.
    void updateDistOnPush(int w, const std::vector<int>& universe,
                          std::vector<std::pair<int,int>>& dist_undo) {
        bfs_touched.clear();
        bfs_queue.clear();
        bfs_d[w] = 0;
        bfs_touched.push_back(w);
        bfs_queue.push_back(w);
        for (size_t qi = 0; qi < bfs_queue.size(); ++qi) {
            int x = bfs_queue[qi];
            int dx = bfs_d[x];
            for (int y : G.adj[x]) {
                if (in_scope[y] && bfs_d[y] == INT_MAX) {
                    bfs_d[y] = dx + 1;
                    bfs_touched.push_back(y);
                    bfs_queue.push_back(y);
                }
            }
        }
        for (int u : universe) {
            if (!in_scope[u] || in_S[u]) continue;
            int sp = bfs_d[u];                         // INT_MAX if unreachable in G[H]
            if (sp > dist_S[u]) { dist_undo.push_back({u, dist_S[u]}); dist_S[u] = sp; }
        }
        for (int x : bfs_touched) bfs_d[x] = INT_MAX;  // reset scratch
    }

    // Undo dist_S raises recorded at index >= `from` (LIFO).
    void restoreDist(std::vector<std::pair<int,int>>& dist_undo, int from) {
        for (int j = (int)dist_undo.size() - 1; j >= from; --j)
            dist_S[dist_undo[j].first] = dist_undo[j].second;
        dist_undo.resize(from);
    }

    // refine: per-node pruning, from the paper's two rules plus the §5.3
    // pivot feasibility check —
    //   Rule 1 (degree): drop candidates with d_a < θ and cascade; prune the node
    //                    if any S vertex has d_a < θ;
    //   Rule 2 (diameter): drop candidates too far from / disconnected from S;
    //   size bound: prune the node if |S| exceeds its size bound;
    //   §5.3 saturation: prune the node if some S vertex already has more
    //                    non-neighbours in S than any feasible flexi-clique tolerates.
    // Mutates Cr/Cun (drops pruned candidates); records removals for undo.
    // Returns false to prune the whole node.
    bool refine(const std::vector<int>& S,
                std::vector<int>& Cr, std::vector<int>& Cun,
                std::vector<std::pair<int,int>>& deg_undo,
                std::vector<int>& removed) {
        // th is used by Rule 1 (degree) and Rule 2 (diameter); only unused when
        // BOTH are ablated.
#if !defined(NO_RULE1) || !defined(NO_DIAMETER)
        int th = theta();
#endif

        // Initial removals: Rule 1 (low adjusted degree) + Rule 2 (diameter).
        // Rule 2 (only when S≠∅) discards a candidate u that would force an
        // infeasible subgraph given its distance du = dist_S[u] from S in G[H],
        // via the two conditions of the diameter rule:
        //   Rule 2(i):  n(θ,du) > |H|              (too many nodes needed to reach u)
        //   Rule 2(ii): ⌊n(θ,du)^τ⌋ > min_{v∈S} d_a(v)  (would over-constrain some S vertex)
        // du = INT_MAX (u disconnected from S in G[H]) is always pruned — this is
        // what lets the rule delete unreachable candidates. Both conditions use the
        // node's incoming |H| and min-degree, evaluated before any cascade.
#ifndef NO_DIAMETER
        const bool Snon = !S.empty();
        const int  scope0 = (int)S.size() + (int)Cr.size() + (int)Cun.size();
        int minAdS = INT_MAX;
        if (Snon) for (int v : S) minAdS = std::min(minAdS, local_adj_deg[v]);
#endif

        std::vector<int> q;
        auto consider = [&](int u) {
            if (!in_scope[u]) return;
#ifdef NO_RULE1
            bool rm = false;                                       // Rule 1 ablated
#else
            bool rm = local_adj_deg[u] < th;                       // Rule 1
#endif
#ifndef NO_DIAMETER
            if (!rm && Snon) {                                     // Rule 2 (diameter)
                int du = dist_S[u];
                if (du == INT_MAX) rm = true;                     // disconnected from S
                else {
                    long long nf = nFuncLL(th, du);
                    if (nf > scope0) rm = true;                                       // Rule 2(i)
                    // Rule 2(ii): ⌊nf^τ⌋ > minAdS  ⟺  nf^τ ≥ minAdS+1. Test the latter
                    // form with a tolerance so double rounding of pow() across an
                    // integer boundary (e.g. 4^0.79248 → 3.0 instead of 2.999…) can
                    // only ever UNDER-prune, never falsely discard a valid extension.
                    else if (std::pow((double)nf, tau) >= (double)minAdS + 1.0 + 1e-6)
                        rm = true;
                }
            }
#endif
            if (rm) q.push_back(u);
        };
        for (int u : Cr)  consider(u);
        for (int u : Cun) consider(u);

        // Cascade (Rule 1 propagation): removing a candidate may drop neighbours
        // below θ, which cascades further.
        for (size_t qi = 0; qi < q.size(); ++qi) {
            int v = q[qi];
            if (!in_scope[v]) continue;          // already removed
            in_scope[v] = false;
            removed.push_back(v);
            for (int u : G.adj[v]) {
                if (!in_scope[u]) continue;
                deg_undo.push_back({u, local_adj_deg[u]});
                --local_adj_deg[u];
#ifndef NO_RULE1
                if (local_adj_deg[u] < th && !in_S[u]) q.push_back(u);   // Rule 1 cascade
#endif
            }
        }
        Cr.erase(std::remove_if(Cr.begin(), Cr.end(),
                                [&](int u){ return !in_scope[u]; }), Cr.end());
        Cun.erase(std::remove_if(Cun.begin(), Cun.end(),
                                 [&](int u){ return !in_scope[u]; }), Cun.end());

        // Rule 1 on S.
#ifndef NO_RULE1
        for (int v : S)
            if (local_adj_deg[v] < th) return false;
#endif

        // Rule 2: size upper bound from min adjusted degree.
        // F_max = ⌊(min_{v∈S} d_a(v) + 1)^{1/τ}⌋ is the largest size any flexi
        // containing S can reach. Prune if |S| already exceeds it, OR |S| == F_max
        // but S is not yet a valid flexi (it can no longer grow, so it never will).
#ifndef NO_RULE2
        if (!S.empty()) {
            int md = INT_MAX;
            for (int v : S) md = std::min(md, local_adj_deg[v]);
            int fmax = floorRoot(md + 1);
            if ((int)S.size() > fmax) return false;
            if ((int)S.size() == fmax && !feasibleS(S)) return false;
        }
#endif

        // §5.3 saturation check (the width-based pivot's feasibility test, not a
        // separate rule): a node is infeasible once some S vertex has more
        // non-neighbours inside S than any feasible flexi-clique can tolerate.
        int scope_size = (int)S.size() + (int)Cr.size() + (int)Cun.size();
        int s_size = (int)S.size();
        for (int v : S) {
            if (smaxOf(v, scope_size) < s_size) return false;
            if (nonDegInS(v, s_size) > dbarMax(v, scope_size)) return false;
        }
        return true;
    }

    // Classify the node's candidate universe into the child's (Cr, Cun) by
    // current state: in-S / out-of-scope are skipped; deg_in_S>0 ⇒ reachable.
    void buildChildSets(const std::vector<int>& cand_uni,
                        std::vector<int>& childCr, std::vector<int>& childCun) {
        childCr.clear();
        childCun.clear();
        for (int u : cand_uni) {
            if (in_S[u] || !in_scope[u]) continue;
            if (deg_in_S[u] > 0) childCr.push_back(u);
            else                 childCun.push_back(u);
        }
    }

    // Candidate ordering for CP-branching = the width-based pivot (§5.3, default ON).
    //
    // Pick the least-slack node p = argmin_{v∈S∪Cr} margin(v), where
    //   margin(v) = dbar_max(v) − dbar_S(v)   (non-neighbour slack until v breaks).
    // Branch ONLY on p (if p∈Cr) plus the first t = min(margin(p), |N̄(p)|)
    // non-neighbours of p in Cr (ascending margin); every other candidate —
    // notably p's neighbours — is NOT branched individually but left in Cr and
    // deferred to the include-all child B_k and deeper recursion. By the Width
    // Bound lemma this caps the exclusion children at ≈margin(p)+1: front-loading
    // p's non-neighbours drives dbar_S(p) to dbar_max(p) so the saturation prune
    // fires, and the deferred candidates cannot extend S without over-saturating
    // p, so completeness (and the maximum) is preserved.
    //
    // NO_PIVOT (ablation) falls back to branching on the full Cr, degree-ordered.
    std::vector<int> orderCr(const std::vector<int>& S,
                             const std::vector<int>& Cr,
                             const std::vector<int>& Cun, int scope_size) {
        (void)Cun;
#ifndef NO_PIVOT
        auto margin = [&](int v) {
            return dbarMax(v, scope_size) - nonDegInS(v, (int)S.size());
        };
        int p = -1, pm = INT_MAX;
        for (int v : S)  { int m = margin(v); if (m < pm) { pm = m; p = v; } }
        for (int v : Cr) { int m = margin(v); if (m < pm) { pm = m; p = v; } }

        // N̄(p) = non-neighbours of p inside Cr, ascending margin.
        std::vector<std::pair<int,int>> non_nbr;
        for (int u : Cr)
            if (u != p && !G.hasEdge(p, u)) non_nbr.push_back({margin(u), u});
        std::sort(non_nbr.begin(), non_nbr.end());

        int t = std::min(std::max(pm, 0), (int)non_nbr.size());  // min(margin(p), |N̄(p)|)
        std::vector<int> order;
        order.reserve(t + 1);
        if (!in_S[p]) order.push_back(p);          // p ∈ Cr → branch on p itself first
        for (int i = 0; i < t; ++i) order.push_back(non_nbr[i].second);

        // Empty ordering (t==0 and p∈S): fall back to the full Cr so the node can
        // still make progress.
        if (order.empty()) {
            order = Cr;
            std::sort(order.begin(), order.end(),
                      [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
        }
        return order;
#else
        // Pivot ablation: branch on the FULL Cr (no width control), ordered by:
        //   default          -> adjusted degree ascending
        //   HIGH_DEGREE_ORDER -> adjusted degree descending
        //   RANDOM_ORDER      -> reproducible random shuffle
        (void)S; (void)scope_size;
        std::vector<int> order = Cr;
#if defined(RANDOM_ORDER)
        std::shuffle(order.begin(), order.end(), rng);
#elif defined(HIGH_DEGREE_ORDER)
        std::sort(order.begin(), order.end(),
                  [&](int a, int b){ return local_adj_deg[a] > local_adj_deg[b]; });
#else
        std::sort(order.begin(), order.end(),
                  [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
#endif
        return order;
#endif
    }

    void recurse(std::vector<int>& S, std::vector<int> Cr, std::vector<int> Cun) {
        ++branches;
#ifdef TRACE_HIST
        ::g_hist[(int)std::min<size_t>(S.size(), 255)]++;
#endif

        // Incumbent update (before refine, as in the oracle).
        if ((int)S.size() > best_size && feasibleS(S)) {
            best = S;
            best_size = (int)S.size();
        }

        std::vector<std::pair<int,int>> ref_deg_undo;
        std::vector<int> ref_removed;
        if (!refine(S, Cr, Cun, ref_deg_undo, ref_removed)) {
            undoScope(ref_removed, ref_deg_undo);
            return;
        }

        // Sound size bound: the whole reachable universe can't beat the
        // incumbent, so neither can any descendant. (max achievable ≤ |R|.)
        if ((int)S.size() + (int)Cr.size() + (int)Cun.size() <= best_size) {
            undoScope(ref_removed, ref_deg_undo);
            return;
        }

        // Seed selection at S = ∅ (§2.5): each Cun vertex seeds one component,
        // earlier seeds are excluded so every connected subgraph is enumerated
        // exactly once by its lowest-ranked vertex.
        if (S.empty()) {
            std::vector<int> seeds = Cun;       // snapshot of refined Cun
            // The root component enumeration dominates the search, so the candidate
            // ordering ablation (default OFF) must order the seeds too, consistently
            // with orderCr(). Every connected subgraph is still enumerated exactly
            // once (by its first-processed seed), so the order affects efficiency,
            // not correctness. Default build = adjusted-degree ascending (unchanged).
#if defined(RANDOM_ORDER)
            std::shuffle(seeds.begin(), seeds.end(), rng);
#elif defined(HIGH_DEGREE_ORDER)
            std::sort(seeds.begin(), seeds.end(),
                      [&](int a, int b){ return local_adj_deg[a] > local_adj_deg[b]; });
#else
            std::sort(seeds.begin(), seeds.end(),
                      [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
#endif
            std::vector<std::pair<int,int>> seed_deg_undo;
            std::vector<int> seed_removed;
            std::vector<int> childCr, childCun;
            std::vector<std::pair<int,int>> seed_dist_undo;
            for (int w : seeds) {
                if (!in_scope[w]) continue;
                pushToS(w, S);
                updateDistOnPush(w, Cun, seed_dist_undo);
                buildChildSets(Cun, childCr, childCun);
                recurse(S, childCr, childCun);
                restoreDist(seed_dist_undo, 0);
                popFromS(w, S);
                in_D[w] = true;
                scopeRemove(w, seed_deg_undo, seed_removed);
            }
            undoScope(seed_removed, seed_deg_undo, /*reset_in_D=*/true);
            undoScope(ref_removed, ref_deg_undo);
            return;
        }

        if (Cr.empty()) {                       // leaf
            undoScope(ref_removed, ref_deg_undo);
            return;
        }

        int scope_size = (int)S.size() + (int)Cr.size() + (int)Cun.size();
        std::vector<int> order = orderCr(S, Cr, Cun, scope_size);

        std::vector<int> cand_uni;
        cand_uni.reserve(Cr.size() + Cun.size());
        cand_uni.insert(cand_uni.end(), Cr.begin(), Cr.end());
        cand_uni.insert(cand_uni.end(), Cun.begin(), Cun.end());

        // Sibling-monotonicity break (§5.2): across siblings S_0⊆S_1⊆…⊆S_k, so
        // the partial set only grows. If the current S (= S_i) is already
        // size-infeasible (Rule 2) or saturated, then every later sibling and the
        // include-all child B_k contains it and is infeasible too — stop. Evaluated
        // in the node's (un-excluded) scope, which is stable under prefix
        // accumulation and conservative w.r.t. each child's smaller scope, so the
        // break is sound. This is what bounds cumulative branching width: without
        // it the size bound never bites (each child excludes only one candidate).
        auto siblingInfeasible = [&]() -> bool {
            int s_size = (int)S.size();
#ifndef NO_RULE2
            int md = INT_MAX;
            for (int v : S) md = std::min(md, local_adj_deg[v]);
            int fmax = floorRoot(md + 1);
            if (s_size > fmax) return true;
            if (s_size == fmax && !feasibleS(S)) return true;
#endif
            for (int v : S) {
                if (smaxOf(v, scope_size) < s_size) return true;
                if (nonDegInS(v, s_size) > dbarMax(v, scope_size)) return true;
            }
            return false;
        };

        std::vector<int> childCr, childCun;
        std::vector<std::pair<int,int>> dist_undo;   // dist_S raises across the prefix
        std::vector<int> dist_marks;                 // dist_undo size before each push
        int pushed = 0;
        bool broke = false;
        for (size_t idx = 0; idx < order.size(); ++idx) {
            if (siblingInfeasible()) { broke = true; break; }
            int vi = order[idx];                // S == nodeS ∪ {order[0..idx-1]}
            // Child B_i: exclude v_i only.
            std::vector<std::pair<int,int>> ex_deg_undo;
            std::vector<int> ex_removed;
            in_D[vi] = true;
            scopeRemove(vi, ex_deg_undo, ex_removed);
            buildChildSets(cand_uni, childCr, childCun);
            recurse(S, childCr, childCun);
            undoScope(ex_removed, ex_deg_undo, /*reset_in_D=*/true);
            // Accumulate v_i into S for later siblings; extend dist_S by a BFS from v_i.
            pushToS(vi, S);
            dist_marks.push_back((int)dist_undo.size());
            updateDistOnPush(vi, cand_uni, dist_undo);
            ++pushed;
        }
        // Include-all child B_k: S == nodeS ∪ order, no exclusion (skipped if the
        // sibling break already proved S_k infeasible).
        if (!broke) {
            buildChildSets(cand_uni, childCr, childCun);
            recurse(S, childCr, childCun);
        }

#ifdef TRACE_WIDTH
        // Width of this CP-branching node = explored exclusion children (`pushed`)
        // plus the include-all child B_k iff the sibling break did not fire.
        {
            int wdt = pushed + (broke ? 0 : 1);
            ::g_width_hist[std::min(wdt, 255)]++;
        }
#endif

        for (int i = pushed - 1; i >= 0; --i) {
            restoreDist(dist_undo, dist_marks[i]);
            popFromS(order[i], S);
        }
        undoScope(ref_removed, ref_deg_undo);
    }
};

// Global degree prune (paper Rule 1, applied globally) on the active mask.
#ifndef NO_RULE1
void globalDegreePrune(const Graph& G, std::vector<bool>& active,
                       std::vector<int>& adeg, int th) {
    std::queue<int> q;
    for (int v = 0; v < G.n; ++v)
        if (active[v] && adeg[v] < th) q.push(v);
    while (!q.empty()) {
        int v = q.front(); q.pop();
        if (!active[v]) continue;
        active[v] = false;
        for (int u : G.adj[v])
            if (active[u] && --adeg[u] < th) q.push(u);
    }
}
#endif

} // namespace

AlgoResult runEBACumulative(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    if (G.n == 0) {                       // empty graph: empty clique
        AlgoResult res;
        res.time_sec = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();
        return res;
    }

    CumState st(G, tau);
    std::vector<int> init_sol;
#ifndef NO_WARMSTART
    // Heuristic initialisation (identical to runEBA): best of FPA / NPA.
    // NO_WARMSTART omits it so the branching comparison isolates the
    // branching scheme from the heuristic lower bound.
    AlgoResult fpa_res = multiSeedGreedyPlusPlus(G, tau);
    AlgoResult npa_res = runNPA(G, tau);
    init_sol = fpa_res.nodes;
    if (npa_res.nodes.size() > init_sol.size()) init_sol = npa_res.nodes;
    if (!init_sol.empty() && isFlexiClique(G, init_sol, tau)) {
        st.best = init_sol;
        st.best_size = (int)init_sol.size();
    }
#endif

    // Global degree prune (paper Rule 1, applied globally) to define V'.
    std::vector<bool> active(G.n, true);
    std::vector<int> adeg(G.n, 0);
    for (int v = 0; v < G.n; ++v) adeg[v] = (int)G.adj[v].size();
#ifndef NO_RULE1
    globalDegreePrune(G, active, adeg, st.theta());
#endif

    // Modified Greedy++ on the pruned subgraph may improve the incumbent.
#ifndef NO_WARMSTART
    {
        std::vector<int> mg = modifiedGreedyPlusPlus(G, tau, &active);
        if ((int)mg.size() > st.best_size && isFlexiClique(G, mg, tau)) {
            st.best = mg;
            st.best_size = (int)mg.size();
            for (int v = 0; v < G.n; ++v) {
                if (!active[v]) continue;
                int d = 0;
                for (int u : G.adj[v]) if (active[u]) ++d;
                adeg[v] = d;
            }
#ifndef NO_RULE1
            globalDegreePrune(G, active, adeg, st.theta());
#endif
        }
    }
#endif  // NO_WARMSTART

    std::vector<int> V_prime;
    for (int v = 0; v < G.n; ++v) if (active[v]) V_prime.push_back(v);

    if (V_prime.empty()) {
        auto t1 = std::chrono::high_resolution_clock::now();
        AlgoResult res;
        res.nodes = st.best.empty() ? init_sol : st.best;
        res.time_sec = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    for (int v : V_prime) {
        st.in_scope[v] = true;
        int d = 0;
        for (int u : G.adj[v]) if (active[u]) ++d;
        st.local_adj_deg[v] = d;
    }

    std::vector<int> S, Cr, Cun = V_prime;
    st.recurse(S, Cr, Cun);

    auto t1 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    res.nodes = st.best.empty() ? init_sol : st.best;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    res.branches = st.branches;
#ifdef TRACE_HIST
    fprintf(stderr, "[ebac depth-hist |S|:count]\n");
    for (int i = 0; i < 256; ++i) if (g_hist[i]) fprintf(stderr, "  %d: %lld\n", i, g_hist[i]);
#endif
#ifdef TRACE_WIDTH
    fprintf(stderr, "[ebac width-hist w:count]\n");
    for (int i = 0; i < 256; ++i) if (g_width_hist[i]) fprintf(stderr, "  %d: %lld\n", i, g_width_hist[i]);
#endif
    return res;
}
