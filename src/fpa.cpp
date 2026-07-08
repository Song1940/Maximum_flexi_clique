#include "fpa.h"
#include "dyn_conn.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <queue>
#include <vector>

// fpa.cpp — Flexi-Prune Algorithm (FPA, Algorithm 1).

// -----------------------------------------------------------
// runFPA — Flexi-Prune Algorithm (Algorithm 1).
//   Step 1: select a dense connected seed from the k-core hierarchy.
//   Step 2: peel the lowest-degree removable node, using a fully dynamic
//           connectivity structure to keep the subgraph connected, until the
//           remainder is a feasible flexi.
// -----------------------------------------------------------
AlgoResult runFPA(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // ── Step 1 ───────────────────────────────────────────────
    std::vector<int> core = G.coreNumbers();
    int k_star = *std::max_element(core.begin(), core.end());

    std::vector<bool> seed_mask(G.n, false);
    int best_k = k_star;

    for (int k = 2; k <= k_star; ++k) {
        std::vector<bool> mask(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= k) mask[v] = true;
        std::vector<bool> lcc = G.largestCC(mask);
        int sz = 0; for (bool b : lcc) sz += b;
        if (floorPow(sz, tau) <= k) { best_k = k; seed_mask = lcc; break; }
    }

    // Scope = (best_k − 1)-core component containing the seed.
    std::vector<bool> scope(G.n, false);
    {
        std::vector<bool> lower(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= best_k - 1) lower[v] = true;

        int seed_node = -1;
        for (int v = 0; v < G.n; ++v) if (seed_mask[v]) { seed_node = v; break; }

        if (seed_node >= 0) {
            std::vector<bool> vis(G.n, false);
            std::queue<int> q;
            q.push(seed_node); vis[seed_node] = true;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                scope[u] = true;
                for (int w : G.adj[u])
                    if (lower[w] && !vis[w]) { vis[w] = true; q.push(w); }
            }
        } else {
            scope = lower;
        }
    }

    // ── Connectivity-preserving peeling with fully dynamic connectivity ──────
    // (Algorithm 1, lines 7-16). The current subgraph is maintained in the
    // dynamic connectivity structure of Holm et al.; a candidate is removable
    // iff deleting its incident edges raises the component count by exactly one
    // (it only isolates itself and does not split the subgraph). We peel the
    // lowest-degree removable node until the remainder is a feasible flexi.
    std::vector<char> inC(G.n, 0);
    std::vector<int>  node_deg(G.n, 0);
    int hsize = 0, max_deg = 0;
    for (int v = 0; v < G.n; ++v) if (scope[v]) { inC[v] = 1; ++hsize; }
    for (int v = 0; v < G.n; ++v)
        if (inC[v]) {
            int d = 0;
            for (int u : G.adj[v]) if (inC[u]) ++d;
            node_deg[v] = d;
            max_deg = std::max(max_deg, d);
        }

    FullyDynamicConnectivity fdc(G.n);
    for (int v = 0; v < G.n; ++v)
        if (inC[v])
            for (int u : G.adj[v])
                if (inC[u] && v < u) fdc.insert_edge(v, u);

    // Incumbent: the k-core seed, used only if it is itself a feasible flexi.
    std::vector<int> seed_nodes;
    for (int v = 0; v < G.n; ++v) if (seed_mask[v]) seed_nodes.push_back(v);
    std::vector<int> best_result;
    int best_size = 0;
    if (isFlexiClique(G, seed_nodes, tau)) { best_result = seed_nodes; best_size = (int)seed_nodes.size(); }

    std::vector<std::vector<int>> buckets(max_deg + 2);
    for (int v = 0; v < G.n; ++v) if (inC[v]) buckets[node_deg[v]].push_back(v);

    auto min_deg_in_C = [&]() -> int {
        for (int d = 0; d <= max_deg; ++d)
            for (int v : buckets[d])
                if (inC[v] && node_deg[v] == d) return d;
        return 0;
    };

    std::vector<char> checked(G.n, 0);
    std::vector<int>  checked_list;
    auto clear_checked = [&]() { for (int v : checked_list) checked[v] = 0; checked_list.clear(); };

    int current_degree = min_deg_in_C();

    while (hsize > best_size) {
        if (min_deg_in_C() >= floorPow(hsize, tau)) {         // current subgraph is a flexi
            best_result.clear();
            for (int v = 0; v < G.n; ++v) if (inC[v]) best_result.push_back(v);
            best_size = hsize;
            break;
        }

        auto bucket_ready = [&](int d) -> bool {
            if (d > max_deg) return false;
            for (int v : buckets[d])
                if (inC[v] && node_deg[v] == d && !checked[v]) return true;
            return false;
        };
        while (current_degree <= max_deg && !bucket_ready(current_degree)) {
            ++current_degree;
            clear_checked();
        }
        if (current_degree > max_deg) break;

        int candidate = -1;
        for (int v : buckets[current_degree])
            if (inC[v] && node_deg[v] == current_degree && !checked[v]) { candidate = v; break; }
        if (candidate < 0) { ++current_degree; clear_checked(); continue; }

        checked[candidate] = 1; checked_list.push_back(candidate);

        long long old_nc = fdc.get_component_count();
        for (int nb : G.adj[candidate]) if (inC[nb]) fdc.delete_edge(candidate, nb);
        long long diff = fdc.get_component_count() - old_nc;

        if (diff == 1) {                                      // non-cut node: remove it
            inC[candidate] = 0;
            --hsize;
            for (int nb : G.adj[candidate])
                if (inC[nb]) { int nd = --node_deg[nb]; if (nd >= 0) buckets[nd].push_back(nb); }
            node_deg[candidate] = 0;
            clear_checked();
            current_degree = min_deg_in_C();
        } else {                                              // cut vertex: roll back
            for (int nb : G.adj[candidate]) if (inC[nb]) fdc.insert_edge(candidate, nb);
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    AlgoResult res;
    res.nodes    = best_result;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    return res;
}
