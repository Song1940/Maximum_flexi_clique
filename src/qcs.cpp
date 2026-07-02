// qcs.cpp — Query-centric Community Search via maximum degree-constrained subgraph.
//
// For a query node q and a "model" (flexi / quasi / kplex), find the MAXIMUM
// vertex set K with q in K such that every vertex of K has internal degree
// >= f(|K|), where the model fixes f(n):
//
//   flexi(tau)  : f(n) = floor(n^tau)              (connectivity REQUIRED)
//   quasi(gamma): f(n) = ceil(gamma * (n - 1))     (connectivity not required)
//   kplex(k)    : f(n) = max(0, n - k)             (connectivity not required)
//
// All three are "min-degree >= f(n)" relaxations with f non-decreasing in n,
// so a single branch-and-bound solver handles them by swapping f.  This makes
// the cross-model comparison maximally fair (identical code path).
//
// The search is restricted to a local region around q (BFS up to --hops,
// capped at --cap closest vertices) and is ANYTIME: a per-query wall-clock
// limit (--timeout) returns the best valid subgraph found so far.  The
// singleton {q} is always valid (f(1)=0), so a non-empty result is guaranteed.
//
// Batch mode: load the graph once, answer many queries from --queries.
//   query line:  <seed_orig> <model> <param>
//   output line: seed=<orig> model=<m> param=<p> size=<s> time=<t> valid=<0|1> nodes=<orig,csv>
//
// Node IDs in BOTH the query file and the output are ORIGINAL ids (as in the
// edge-list file); the internal re-indexing is hidden.

#include "graph.h"
#include "eba.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using Clock = std::chrono::steady_clock;

// ─── Graph loader that PRESERVES the original node ids ──────────────────────
struct LabeledGraph {
    Graph G;
    std::vector<int>           orig;        // internal -> original id
    std::unordered_map<int,int> to_internal; // original id -> internal
};

static LabeledGraph loadLabeled(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot open: " + path);

    std::unordered_map<int,int> id_map;
    std::vector<int>            orig;
    std::vector<std::pair<int,int>> edges;
    int next_id = 0;

    auto getID = [&](int raw) -> int {
        auto [it, inserted] = id_map.emplace(raw, next_id);
        if (inserted) { ++next_id; orig.push_back(raw); }
        return it->second;
    };

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        long long u, v;
        if (!(ss >> u >> v) || u == v) continue;
        edges.push_back({getID((int)u), getID((int)v)});
    }

    Graph G(next_id);
    for (auto [u, v] : edges) G.addEdge(u, v);
    G.finalize();

    LabeledGraph LG;
    LG.G = std::move(G);
    LG.orig = std::move(orig);
    LG.to_internal = std::move(id_map);
    std::cerr << "Loaded graph: n=" << LG.G.n << " m=" << LG.G.m << "\n";
    return LG;
}

// ─── Model degree-threshold function f(n) ──────────────────────────────────
enum class Model { Flexi, Quasi, KPlex };

struct ModelSpec {
    Model  model;
    double param;          // tau, gamma, or k
    bool   need_connected; // flexi only
    // Required minimum internal degree for a solution of size n.
    int f(int n) const {
        if (n <= 1) return 0;
        switch (model) {
            case Model::Flexi:
                // Route through the shared floorPow so the case study's flexi
                // threshold is bit-identical to the paper's solver definition.
                return floorPow(n, param);
            case Model::Quasi:
                // Subtract a tiny epsilon before ceil to absorb float rounding
                // (e.g. 0.7*10 == 7.0000001 must not ceil up to 8).
                return (int)std::ceil(param * (double)(n - 1) - 1e-9);
            case Model::KPlex: {
                int k = (int)std::lround(param);
                return std::max(0, n - k);
            }
        }
        return 0;
    }
};

// ─── Local solver: works on a compact region graph ──────────────────────────
// Region vertices are relabeled to [0, R).  radj[i] = sorted local neighbors.
//
// State is held IN PLACE (no per-node vector copies) and mutated through an
// undo stack, so backtracking is O(changes) not O(R):
//   status[v] : 2 = inset (forced in), 1 = cand (undecided), 0 = removed
//   potdeg[v] : # neighbours of v that are NOT removed (inset∪cand)
//               = optimistic upper bound on v's achievable internal degree.
// Removing a vertex decrements potdeg of its live neighbours; undo reverses it.
struct LocalSolver {
    const std::vector<std::vector<int>>& radj;
    int R;
    int seed;                 // local index of the query vertex
    ModelSpec spec;
    Clock::time_point deadline;
    bool timed_out = false;

    std::vector<int> best;    // local indices of best solution found
    long long node_budget = 0;

    std::vector<char> status; // 2 inset, 1 cand, 0 removed
    std::vector<int>  potdeg;  // live-neighbour count
    int n_in = 0, n_cand = 0;

    struct Op { int type; int v; };   // type 0 = was-removed (cand→0), 1 = was-included (cand→2)
    std::vector<Op> undo;
    std::vector<char> seen;   // reusable BFS buffer

    LocalSolver(const std::vector<std::vector<int>>& radj_, int seed_,
                ModelSpec spec_, Clock::time_point dl)
        : radj(radj_), R((int)radj_.size()), seed(seed_), spec(spec_), deadline(dl),
          status(R, 1), potdeg(R, 0), seen(R, 0) {}

    bool outOfTime() {
        if ((++node_budget & 1023) == 0 && Clock::now() >= deadline) timed_out = true;
        return timed_out;
    }

    // cand → removed.
    void doRemove(int v) {
        status[v] = 0; --n_cand;
        for (int u : radj[v]) if (status[u]) --potdeg[u];
        undo.push_back({0, v});
    }
    // cand → inset (potdeg unaffected: v stays live).
    void doInclude(int v) {
        status[v] = 2; --n_cand; ++n_in;
        undo.push_back({1, v});
    }
    void rewind(size_t cp) {
        while (undo.size() > cp) {
            Op op = undo.back(); undo.pop_back();
            if (op.type == 0) {               // undo a removal
                status[op.v] = 1; ++n_cand;
                for (int u : radj[op.v]) if (status[u]) ++potdeg[u];
            } else {                          // undo an inclusion
                status[op.v] = 1; --n_in; ++n_cand;
            }
        }
    }

    // Remove cand vertices not reachable from seed through (inset∪cand).
    bool connectivityPrune() {
        std::fill(seen.begin(), seen.end(), 0);
        std::queue<int> q;
        q.push(seed); seen[seed] = 1;         // seed is always inset
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int w : radj[u])
                if (status[w] && !seen[w]) { seen[w] = 1; q.push(w); }
        }
        bool changed = false;
        for (int v = 0; v < R; ++v)
            if (status[v] == 1 && !seen[v]) { doRemove(v); changed = true; }
        return changed;
    }

    // Reduce to fixpoint.  Returns false if this branch is provably dead.
    bool reduce() {
        for (;;) {
            int target = std::max(n_in, (int)best.size() + 1);
            int thr = spec.f(target);
            bool changed = false;
            for (int v = 0; v < R; ++v) {
                if (status[v] == 2) {                 // forced vertex
                    if (potdeg[v] < thr) return false;
                } else if (status[v] == 1) {          // candidate
                    if (potdeg[v] < thr) { doRemove(v); changed = true; }
                }
            }
            if (spec.need_connected && connectivityPrune()) changed = true;
            if (n_in + n_cand <= (int)best.size()) return false;
            if (!changed) break;
        }
        return true;
    }

    // Upper bound on achievable solution size from the current state.
    // For a size-s solution every member needs internal degree ≥ f(s) ≤ potdeg,
    // so we need ≥ s live vertices with potdeg ≥ f(s), and every inset vertex
    // must reach f(s).  Return the largest feasible s (a necessary condition).
    int sizeUpperBound() {
        int N = n_in + n_cand;
        // suffix[d] = # live vertices with potdeg ≥ d  (potdeg ≤ R).
        std::vector<int> cnt(R + 2, 0);
        int min_in = INT_MAX;
        for (int v = 0; v < R; ++v) {
            if (!status[v]) continue;
            ++cnt[potdeg[v]];
            if (status[v] == 2) min_in = std::min(min_in, potdeg[v]);
        }
        for (int d = R; d >= 0; --d) cnt[d] += cnt[d + 1];   // now cnt[d] = #(potdeg ≥ d)
        if (min_in == INT_MAX) min_in = R;                   // no inset (shouldn't happen)
        for (int s = N; s >= n_in; --s) {
            int need = spec.f(s);
            if (need > min_in) continue;                     // some inset vertex can't reach
            if (cnt[std::min(need, R + 1)] >= s) return s;
        }
        return n_in;
    }

    // Internal-degree validity of the current inset.
    bool isValidInset() {
        if (n_in <= 1) return true;
        int thr = spec.f(n_in);
        for (int v = 0; v < R; ++v) {
            if (status[v] != 2) continue;
            int d = 0;
            for (int u : radj[v]) if (status[u] == 2) ++d;
            if (d < thr) return false;
        }
        if (spec.need_connected) {
            std::fill(seen.begin(), seen.end(), 0);
            std::queue<int> q;
            q.push(seed); seen[seed] = 1; int c = 1;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int w : radj[u])
                    if (status[w] == 2 && !seen[w]) { seen[w] = 1; ++c; q.push(w); }
            }
            if (c != n_in) return false;
        }
        return true;
    }

    void branch() {
        if (outOfTime()) return;
        size_t cp = undo.size();

        if (!reduce()) { rewind(cp); return; }

        if (n_in > (int)best.size() && isValidInset()) {
            best.clear();
            for (int v = 0; v < R; ++v) if (status[v] == 2) best.push_back(v);
        }

        if (n_cand == 0) { rewind(cp); return; }
        if (sizeUpperBound() <= (int)best.size()) { rewind(cp); return; }

        // Branch vertex: live candidate with smallest optimistic degree.
        int u = -1, ub = INT_MAX;
        for (int v = 0; v < R; ++v)
            if (status[v] == 1 && potdeg[v] < ub) { ub = potdeg[v]; u = v; }
        if (u == -1) { rewind(cp); return; }

        size_t cpb = undo.size();
        doInclude(u); branch(); rewind(cpb);     // include u
        if (timed_out) { rewind(cp); return; }
        doRemove(u);  branch();                  // exclude u
        rewind(cp);
    }

    std::vector<int> solve() {
        best = { seed };
        for (int v = 0; v < R; ++v) potdeg[v] = (int)radj[v].size();
        status.assign(R, 1);
        status[seed] = 2;
        n_in = 1; n_cand = R - 1;
        branch();
        return best;
    }
};

// ─── Build the local region around a seed and solve ─────────────────────────
struct QueryResult {
    int  size = 0;
    bool valid = false;
    bool timed_out = false;   // true => best-so-far (not proven optimal within region)
    double time_sec = 0;
    std::vector<int> nodes_internal;
};

static QueryResult runQuery(const Graph& G, int seed, ModelSpec spec,
                            int hops, int cap, double timeout_sec,
                            const std::vector<int>* restrict_set = nullptr) {
    auto t0 = Clock::now();

    // The search region is either an explicit vertex set (coverage mode: the
    // ground-truth community, so predictions ⊆ C and precision is fixed at 1)
    // or the BFS ball of radius `hops` around the seed (recovery mode).
    std::vector<int> region;
    std::unordered_map<int,int> local_of;     // global -> local index
    if (restrict_set) {
        // Seed first so it always has a local index, then the rest of the set.
        region.push_back(seed); local_of[seed] = 0;
        for (int v : *restrict_set)
            if (v != seed && !local_of.count(v)) {
                local_of[v] = (int)region.size();
                region.push_back(v);
            }
    } else {
        std::queue<std::pair<int,int>> q;      // (vertex, distance)
        std::unordered_map<int,int> dist;
        q.push({seed, 0}); dist[seed] = 0;
        while (!q.empty() && (int)region.size() < cap) {
            auto [v, d] = q.front(); q.pop();
            local_of[v] = (int)region.size();
            region.push_back(v);
            if (d == hops) continue;
            for (int w : G.adj[v])
                if (!dist.count(w)) { dist[w] = d + 1; q.push({w, d + 1}); }
        }
    }

    const int R = (int)region.size();
    std::vector<std::vector<int>> radj(R);
    for (int i = 0; i < R; ++i) {
        int v = region[i];
        for (int w : G.adj[v]) {
            auto it = local_of.find(w);
            if (it != local_of.end()) radj[i].push_back(it->second);
        }
        std::sort(radj[i].begin(), radj[i].end());
    }

    int seed_local = local_of[seed];

    std::vector<int> sol_local;
    bool tout = false;
    if (spec.model == Model::Flexi) {
        // flexi: use the paper's exact EBA (seeded) — far stronger pruning than
        // the generic B&B on the hard lenient-threshold (low-tau) regime.
        Graph rg(R);
        for (int i = 0; i < R; ++i)
            for (int j : radj[i]) if (j > i) rg.addEdge(i, j);
        rg.finalize();
        AlgoResult er = runEBASeeded(rg, seed_local, spec.param, timeout_sec);
        sol_local = er.nodes;
        tout = er.timed_out;
    } else {
        // quasi / kplex: unified anytime branch-and-bound.
        auto deadline = t0 + std::chrono::duration_cast<Clock::duration>(
                                std::chrono::duration<double>(timeout_sec));
        LocalSolver solver(radj, seed_local, spec, deadline);
        sol_local = solver.solve();
        tout = solver.timed_out;
    }

    QueryResult res;
    res.size = (int)sol_local.size();
    res.timed_out = tout;
    for (int li : sol_local) res.nodes_internal.push_back(region[li]);

    // Independent validity re-check against the ORIGINAL graph (defensive).
    // Uses a hash-set over the (small) result, so cost is O(sum-deg of result)
    // rather than O(G.n) per query — important on million-node graphs.
    res.valid = true;
    {
        int n = res.size;
        if (n > 1) {
            int thr = spec.f(n);
            std::unordered_set<int> mark(res.nodes_internal.begin(),
                                         res.nodes_internal.end());
            for (int v : res.nodes_internal) {
                int d = 0;
                for (int u : G.adj[v]) if (mark.count(u)) ++d;
                if (d < thr) { res.valid = false; break; }
            }
            if (res.valid && spec.need_connected) {
                std::unordered_set<int> seen;
                std::queue<int> bq;
                int s0 = res.nodes_internal[0];
                bq.push(s0); seen.insert(s0);
                while (!bq.empty()) {
                    int u = bq.front(); bq.pop();
                    for (int w : G.adj[u])
                        if (mark.count(w) && !seen.count(w)) { seen.insert(w); bq.push(w); }
                }
                if ((int)seen.size() != n) res.valid = false;
            }
        }
    }
    res.time_sec = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

static bool parseModel(const std::string& s, double param, ModelSpec& out) {
    if (s == "flexi") { out = {Model::Flexi, param, true};  return true; }
    if (s == "quasi") { out = {Model::Quasi, param, false}; return true; }
    if (s == "kplex") { out = {Model::KPlex, param, false}; return true; }
    return false;
}

static void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --file <network.dat> --queries <file> [--out <file>]"
              << " [--hops 2] [--cap 3000] [--timeout 5.0]\n"
              << "  query line: <seed_orig> <flexi|quasi|kplex> <param>\n";
}

int main(int argc, char** argv) {
    std::string file_path, query_path, out_path;
    int hops = 2, cap = 3000;
    double timeout_sec = 5.0;

    for (int i = 1; i < argc; ++i) {
        if      (!std::strcmp(argv[i], "--file")    && i+1 < argc) file_path  = argv[++i];
        else if (!std::strcmp(argv[i], "--queries") && i+1 < argc) query_path = argv[++i];
        else if (!std::strcmp(argv[i], "--out")     && i+1 < argc) out_path   = argv[++i];
        else if (!std::strcmp(argv[i], "--hops")    && i+1 < argc) hops       = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--cap")     && i+1 < argc) cap        = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--timeout") && i+1 < argc) timeout_sec= std::atof(argv[++i]);
    }
    if (file_path.empty() || query_path.empty()) { printUsage(argv[0]); return 1; }
    if (cap  < 1) cap  = 1;     // empty region would UB the solver
    if (hops < 0) hops = 0;
    if (timeout_sec <= 0) timeout_sec = 0.001;

    LabeledGraph LG;
    try { LG = loadLabeled(file_path); }
    catch (const std::exception& e) { std::cerr << "Error: " << e.what() << "\n"; return 1; }

    std::ifstream qin(query_path);
    if (!qin.is_open()) { std::cerr << "Cannot open queries: " << query_path << "\n"; return 1; }

    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (!out_path.empty()) { fout.open(out_path); out = &fout; }

    std::string line;
    int processed = 0, skipped = 0;
    while (std::getline(qin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        long long seed_orig; std::string model_str; double param;
        if (!(ss >> seed_orig >> model_str >> param)) { ++skipped; continue; }

        // Optional 4th token: comma-separated restriction set (original ids).
        // When present the search is confined to these vertices (coverage mode:
        // pass the ground-truth community so predictions ⊆ C ⇒ precision = 1).
        std::vector<int> restrict_internal;
        std::string rtok;
        bool has_restrict = static_cast<bool>(ss >> rtok);

        ModelSpec spec;
        if (!parseModel(model_str, param, spec)) { ++skipped; continue; }

        auto it = LG.to_internal.find((int)seed_orig);
        if (it == LG.to_internal.end()) {
            // Seed not present in the graph (isolated / missing): emit singleton.
            (*out) << "seed=" << seed_orig << " model=" << model_str
                   << " param=" << param << " size=0 time=0 valid=0 timeout=0 nodes=\n";
            ++skipped; continue;
        }
        int seed = it->second;

        if (has_restrict) {
            std::stringstream rs(rtok);
            std::string tok;
            while (std::getline(rs, tok, ',')) {
                if (tok.empty()) continue;
                auto rit = LG.to_internal.find(std::atoi(tok.c_str()));
                if (rit != LG.to_internal.end()) restrict_internal.push_back(rit->second);
            }
        }

        QueryResult res = runQuery(LG.G, seed, spec, hops, cap, timeout_sec,
                                   has_restrict ? &restrict_internal : nullptr);

        (*out) << "seed=" << seed_orig << " model=" << model_str
               << " param=" << param << " size=" << res.size
               << " time=" << res.time_sec << " valid=" << (res.valid ? 1 : 0)
               << " timeout=" << (res.timed_out ? 1 : 0)
               << " nodes=";
        for (size_t i = 0; i < res.nodes_internal.size(); ++i) {
            (*out) << LG.orig[res.nodes_internal[i]];
            if (i + 1 < res.nodes_internal.size()) (*out) << ",";
        }
        (*out) << "\n";
        out->flush();
        ++processed;
    }
    std::cerr << "Processed " << processed << " queries, skipped " << skipped << "\n";
    return 0;
}
