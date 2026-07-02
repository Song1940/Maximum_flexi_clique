#include "eba.h"
#include "fpa.h"
#include "fpa3.h"
#include "graph.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <queue>
#include <vector>

namespace {

// floor(x^{1/tau}), clamped to INT_MAX (small tau -> huge exponent overflows int).
int floorRootB(int x, double tau) {
    if (x <= 0) return 0;
    double r = std::floor(std::pow((double)x, 1.0 / tau));
    if (r >= (double)INT_MAX) return INT_MAX;
    return (int)r;
}

struct BaselineSolver {
    const Graph& G;
    double tau;
    std::vector<int> best;
    int best_size = 0;
    long long branches = 0;
    long long conn_checks = 0;

    BaselineSolver(const Graph& G_, double tau_) : G(G_), tau(tau_) {}

    int theta() const { return computeTheta(best_size, tau); }

    // Degree of v restricted to the marked set.
    int degIn(int v, const std::vector<char>& mark) const {
        int d = 0;
        for (int u : G.adj[v]) if (mark[u]) ++d;
        return d;
    }

    // Is S a valid flexi-clique?  Degree condition + EXPLICIT connectivity check
    // (the overhead that CP avoids).  sm marks S.
    void tryUpdate(const std::vector<int>& S, const std::vector<char>& sm) {
        int k = (int)S.size();
        if (k <= best_size) return;
        int thr = floorPow(k, tau);
        for (int v : S) if (degIn(v, sm) < thr) return;
        ++conn_checks;
        if (k > 1) {
            std::vector<char> seen(G.n, 0);
            std::queue<int> q;
            q.push(S[0]); seen[S[0]] = 1; int c = 1;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int w : G.adj[u])
                    if (sm[w] && !seen[w]) { seen[w] = 1; ++c; q.push(w); }
            }
            if (c != k) return;          // disconnected -> not a valid flexi-clique
        }
        best = S; best_size = k;
    }

    bool prune(const std::vector<int>& S, std::vector<int>& cand,
               std::vector<char>& sm) {
        sm.assign(G.n, 0);
        for (int v : S) sm[v] = 1;
        std::vector<char> scope(G.n, 0);
        for (int v : S)    scope[v] = 1;
        for (int v : cand) scope[v] = 1;

        int th = theta();

        // Rule 1 on S: a forced vertex that can't reach theta kills the node.
        for (int v : S) if (degIn(v, scope) < th) return false;

        if (!S.empty()) {
            int md = INT_MAX;
            for (int v : S) md = std::min(md, degIn(v, scope));
            int fmax = floorRootB(md + 1, tau);
            if ((int)S.size() > fmax) return false;
        }

        // Rule 1 on candidates: drop u with adjusted degree < theta, to fixpoint.
        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<int> keep;
            keep.reserve(cand.size());
            for (int u : cand) {
                if (degIn(u, scope) >= th) keep.push_back(u);
                else { scope[u] = 0; changed = true; }
            }
            cand.swap(keep);
        }

        // Incumbent size bound: no extension here can beat the current best.
        if ((int)S.size() + (int)cand.size() <= best_size) return false;

        // Saturation (non-degree) prune — IDENTICAL to CP's, so the comparison
        // isolates the branching scheme (not the pruning).  For v in S:
        //   smax(v)     = max(1, min(floor((d_a+1)^{1/tau}), scope))   max size w/ v
        //   dbar_max(v) = smax - 1 - floor(smax^tau)        max non-neighbours of v
        //   dbar_S(v)   = (|S|-1) - deg_in_S(v)             current non-neighbours
        // Prune if v can no longer reach |S| (smax < |S|) or already over its
        // non-neighbour budget (dbar_S > dbar_max).
        {
            int scope_size = (int)S.size() + (int)cand.size();
            int ssz = (int)S.size();
            for (int v : S) {
                int da = degIn(v, scope);
                int smax = std::max(1, std::min(floorRootB(da + 1, tau), scope_size));
                int dbar_max = smax - 1 - floorPow(smax, tau);
                int dbar_S = (ssz - 1) - degIn(v, sm);
                if (smax < ssz || dbar_S > dbar_max) return false;
            }
        }
        return true;
    }

    // ── Naive binary include/exclude ────────────────────────────────────────
    void naive(std::vector<int> S, std::vector<int> cand) {
        ++branches;
        std::vector<char> sm;
        if (!prune(S, cand, sm)) return;
        tryUpdate(S, sm);
        if (cand.empty()) return;
        int v = cand.back();
        cand.pop_back();                       // `cand` is now the "rest"
        std::vector<int> Sin = S; Sin.push_back(v);
        naive(std::move(Sin), cand);           // include v
        naive(std::move(S), std::move(cand));  // exclude v
    }

    // ── Set-enumeration (FastQC/IterQC-style) ───────────────────────────────
    void se(std::vector<int> S, std::vector<int> cand) {
        ++branches;
        std::vector<char> sm;
        if (!prune(S, cand, sm)) return;
        tryUpdate(S, sm);
        for (size_t i = 0; i < cand.size(); ++i) {
            std::vector<int> Sin = S; Sin.push_back(cand[i]);
            std::vector<int> later(cand.begin() + i + 1, cand.end());
            se(std::move(Sin), std::move(later));
        }
    }
};

// Global degree prune (paper Rule 1, global form) on a bool active mask: drop
void globalDegPrune(const Graph& G, std::vector<bool>& active, int theta0) {
    std::vector<int> deg(G.n, 0);
    std::queue<int> q;
    for (int v = 0; v < G.n; ++v) {
        if (!active[v]) continue;
        for (int u : G.adj[v]) if (active[u]) ++deg[v];
        if (deg[v] < theta0) q.push(v);
    }
    while (!q.empty()) {
        int v = q.front(); q.pop();
        if (!active[v]) continue;
        active[v] = false;
        for (int u : G.adj[v])
            if (active[u] && --deg[u] < theta0) q.push(u);
    }
}

AlgoResult runBaseline(const Graph& G, double tau, bool set_enum) {
    auto t0 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    if (G.n == 0) return res;

    BaselineSolver st(G, tau);
    std::vector<bool> active(G.n, true);
#ifndef NO_WARMSTART

    AlgoResult fpa = runFPA3(G, tau);
    AlgoResult npa = runNPA(G, tau);
    std::vector<int> init = fpa.nodes;
    if (npa.nodes.size() > init.size()) init = npa.nodes;
    if (!init.empty() && isFlexiClique(G, init, tau)) {
        st.best = init; st.best_size = (int)init.size();
    }
#endif
    globalDegPrune(G, active, st.theta());          // Rule 1 (global)
#ifndef NO_WARMSTART
    {
        std::vector<int> mg = modifiedGreedyPlusPlus(G, tau, &active);
        if ((int)mg.size() > st.best_size && isFlexiClique(G, mg, tau)) {
            st.best = mg; st.best_size = (int)mg.size();
            globalDegPrune(G, active, st.theta());
        }
    }
#endif

    std::vector<int> V_prime;
    for (int v = 0; v < G.n; ++v) if (active[v]) V_prime.push_back(v);
    if (!V_prime.empty()) {
        std::vector<int> S;
        if (set_enum) st.se(S, V_prime);
        else          st.naive(S, V_prime);
    }

    res.nodes = st.best;
    res.branches = st.branches;
    res.depth_prunes = st.conn_checks;   // reused field: # explicit connectivity checks
    res.time_sec = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t0).count();
    return res;
}

} // namespace

AlgoResult runNaiveBranching(const Graph& G, double tau) { return runBaseline(G, tau, false); }
AlgoResult runSEBranching   (const Graph& G, double tau) { return runBaseline(G, tau, true);  }
