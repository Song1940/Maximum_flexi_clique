#include "npa.h"
#include <algorithm>
#include <chrono>

// NPA: peeling-based heuristic baseline for flexi-clique search.
//   1. K-core decomposition to find a seed subgraph.
//   2. Iteratively remove the lowest-degree non-articulation-point node
//      until all nodes satisfy the flexi-clique degree constraint.
AlgoResult runNPA(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Step 1: k-core seed selection.
    std::vector<int> core = G.coreNumbers();
    int k_star = *std::max_element(core.begin(), core.end());

    std::vector<bool> active(G.n, false);
    int best_k = 1;

    for (int k = 2; k <= k_star; ++k) {
        std::vector<bool> mask(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= k) mask[v] = true;

        std::vector<bool> lcc = G.largestCC(mask);
        int sz = 0;
        for (bool b : lcc) sz += b;

        if (floorPow(sz, tau) <= k) {
            best_k = k;
            active = lcc;
            break;
        }
    }

    if (best_k == 1) {
        std::vector<bool> mask(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= k_star) mask[v] = true;
        active = G.largestCC(mask);
    }

    // Expand to the (best_k - 1)-core component that contains the seed.
    {
        std::vector<bool> lower(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= best_k - 1) lower[v] = true;

        std::vector<bool> visited(G.n, false);
        std::queue<int> q;
        for (int v = 0; v < G.n; ++v)
            if (active[v]) { q.push(v); visited[v] = true; break; }

        std::vector<bool> component(G.n, false);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            component[u] = true;
            for (int w : G.adj[u])
                if (lower[w] && !visited[w]) { visited[w] = true; q.push(w); }
        }
        active = component;
    }

    std::vector<int> loc_deg(G.n, 0);
    int active_cnt = 0;
    for (int v = 0; v < G.n; ++v) {
        if (!active[v]) continue;
        ++active_cnt;
        for (int u : G.adj[v])
            if (active[u]) ++loc_deg[v];
    }

    // Step 2: iterative peeling; track the best feasible subgraph seen.
    std::vector<int> best_snapshot;
    int best_snapshot_size = 0;

    auto checkAndSnapshot = [&]() {
        if (active_cnt <= best_snapshot_size) return;
        int thr = floorPow(active_cnt, tau);
        for (int v = 0; v < G.n; ++v)
            if (active[v] && loc_deg[v] < thr) return;
        if (!G.isConnected(active)) return;
        best_snapshot_size = active_cnt;
        best_snapshot.clear();
        for (int v = 0; v < G.n; ++v)
            if (active[v]) best_snapshot.push_back(v);
    };

    checkAndSnapshot();

    bool can_remove = true;
    while (can_remove && active_cnt > 1) {
        can_remove = false;
        int threshold = floorPow(active_cnt, tau);

        // Stop if every node already meets the degree requirement.
        int min_deg = INT_MAX;
        for (int v = 0; v < G.n; ++v)
            if (active[v]) min_deg = std::min(min_deg, loc_deg[v]);
        if (min_deg >= threshold) break;

        // Compute all articulation points once per iteration: O(n+m).
        std::vector<bool> ap = G.articulationPoints(active);

        // Remove the min-degree non-articulation node (no threshold filter).
        int best_u = -1, best_deg = INT_MAX;
        for (int v = 0; v < G.n; ++v) {
            if (!active[v] || ap[v]) continue;
            if (loc_deg[v] < best_deg) { best_deg = loc_deg[v]; best_u = v; }
        }
        if (best_u == -1) break;  // all removable nodes are APs — stuck

        active[best_u] = false;
        --active_cnt;
        for (int w : G.adj[best_u])
            if (active[w]) --loc_deg[w];
        loc_deg[best_u] = 0;
        can_remove = true;
        checkAndSnapshot();
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    AlgoResult res;
    res.nodes    = best_snapshot;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    return res;
}
