#include "fpa.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <queue>
#include <utility>
#include <vector>

// ============================================================
// fpa.cpp — OPTIMISED FPA (lazy-heap MGPP + multi-CC step 3)
//
// Changes vs. the previous version:
//   * modifiedGreedyPlusPlus: O((n+m) log n) instead of O(n^2)
//       - replaces std::unordered_map / std::set with std::vector
//       - replaces the per-peel linear scan with a lazy min-heap
//         (score = load + cur_deg; score is monotone non-increasing)
//       - accepts an optional active_mask so EBA can reuse it on the
//         Rule-6–pruned subgraph without rebuilding a sub-Graph
//   * runFPA step 3: iterates a work-queue over EVERY sufficiently
//     large connected component produced by a batch split, instead
//     of discarding all but the largest.  This is the change that
//     lets FPA find larger solutions on synthetic graphs whose best
//     flexi-clique isn't in the largest post-split CC.
// ============================================================

// -----------------------------------------------------------
// Modified Greedy++ (lazy-heap).
//
// Peels the node with the smallest (load + cur_deg) score until the
// min-score node meets the flexi-clique degree threshold; records the
// largest valid snapshot encountered.  When a threshold-feasible
// snapshot disconnects the graph, each connected component is
// re-queued for independent processing.
// -----------------------------------------------------------
std::vector<int> modifiedGreedyPlusPlus(const Graph& G, double tau,
                                        const std::vector<bool>* active_mask) {
    std::vector<int> best;
    int global_best = 0;

    // Build initial connected components, restricted to active_mask if given.
    std::vector<bool> visited(G.n, false);
    if (active_mask) {
        for (int v = 0; v < G.n; ++v)
            if (!(*active_mask)[v]) visited[v] = true;
    }

    std::vector<std::vector<int>> init_components;
    for (int s = 0; s < G.n; ++s) {
        if (visited[s]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(s); visited[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int w : G.adj[u])
                if (!visited[w]) { visited[w] = true; q.push(w); }
        }
        init_components.push_back(std::move(comp));
    }

    std::sort(init_components.begin(), init_components.end(),
              [](const auto& a, const auto& b){ return a.size() > b.size(); });

    std::queue<std::vector<int>> comp_q;
    for (auto& c : init_components) comp_q.push(std::move(c));

    // Scratch arrays sized G.n, reused across components.
    std::vector<int>    cur_deg(G.n, 0);
    std::vector<double> load(G.n, 0.0);
    std::vector<char>   in_comp(G.n, 0);

    while (!comp_q.empty()) {
        std::vector<int> comp = std::move(comp_q.front());
        comp_q.pop();
        int comp_sz = (int)comp.size();
        if (comp_sz <= global_best) continue;

        // Initialise in_comp / cur_deg / load.
        for (int v : comp) { in_comp[v] = 1; cur_deg[v] = 0; load[v] = 0.0; }
        for (int v : comp)
            for (int u : G.adj[v])
                if (in_comp[u]) ++cur_deg[v];

        int comp_min_deg = INT_MAX;
        for (int v : comp) comp_min_deg = std::min(comp_min_deg, cur_deg[v]);
        int max_possible = (comp_min_deg > 0)
            ? (int)std::floor(std::pow((double)comp_min_deg, 1.0 / tau)) : 0;
        if (max_possible <= global_best) {
            for (int v : comp) in_comp[v] = 0;
            continue;
        }

        // Lazy min-heap.  Because score = load + cur_deg is monotone
        // non-increasing for any fixed v, any popped entry whose score
        // exceeds the current score is stale and can be discarded.
        using Entry = std::pair<double,int>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;
        for (int v : comp) pq.emplace(load[v] + cur_deg[v], v);

        int remaining_cnt = comp_sz;

        while (remaining_cnt > 0) {
            int sz        = remaining_cnt;
            int threshold = floorPow(sz, tau);

            // Extract min-score node (skip stale/removed entries).
            int    min_node = -1;
            while (!pq.empty()) {
                auto top = pq.top();
                double s  = top.first;
                int    v  = top.second;
                if (!in_comp[v]) { pq.pop(); continue; }
                double cur_sc = load[v] + cur_deg[v];
                if (s > cur_sc + 1e-12) { pq.pop(); continue; }  // stale
                min_node = v;
                break;
            }
            if (min_node == -1) break;

            if (cur_deg[min_node] >= threshold) {
                // Snapshot: collect current member list.
                std::vector<int> cur_list;
                cur_list.reserve(sz);
                for (int v : comp) if (in_comp[v]) cur_list.push_back(v);

                // Fast connectivity test via reusable rem_act mask.
                std::vector<bool> rem_act(G.n, false);
                for (int v : cur_list) rem_act[v] = true;

                if (G.isConnected(rem_act)) {
                    if (sz > global_best) {
                        global_best = sz;
                        best        = cur_list;
                    }
                } else {
                    // Split → enqueue every sufficiently-large CC.
                    std::vector<char> cc_vis(G.n, 0);
                    for (int s0 : cur_list) {
                        if (cc_vis[s0]) continue;
                        std::vector<int> sub;
                        std::queue<int> bq;
                        bq.push(s0); cc_vis[s0] = 1;
                        while (!bq.empty()) {
                            int u = bq.front(); bq.pop();
                            sub.push_back(u);
                            for (int w : G.adj[u])
                                if (in_comp[w] && !cc_vis[w]) {
                                    cc_vis[w] = 1; bq.push(w);
                                }
                        }
                        if ((int)sub.size() > global_best)
                            comp_q.push(std::move(sub));
                    }
                }
                break;  // stop peeling this component
            }

            // Peel min_node.
            pq.pop();
            load[min_node] += cur_deg[min_node];
            for (int u : G.adj[min_node]) {
                if (in_comp[u]) {
                    --cur_deg[u];
                    pq.emplace(load[u] + cur_deg[u], u);
                }
            }
            in_comp[min_node] = 0;
            --remaining_cnt;
        }

        for (int v : comp) in_comp[v] = 0;
    }

    return best;
}

// -----------------------------------------------------------
// runFPA
//   Step 1: k-core seed selection.
//   Step 2: Modified Greedy++ on the seed scope.
//   Step 3: connectivity-aware batch peeling across a work-queue
//           of components (NEW: every large CC produced by a split
//           is processed, not just the largest).
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

    // ── Step 2: MGPP seed ────────────────────────────────────
    // Run MGPP both on the scope (fast, aligned with step-3 peeling) and on
    // the full graph (catches solutions whose densest region is not inside
    // the chosen k-core scope — exactly the failure mode we were hitting on
    // synthetic LFR graphs).  Take the larger of the two.
    std::vector<int> greedy_scope = modifiedGreedyPlusPlus(G, tau, &scope);
    std::vector<int> greedy_full  = modifiedGreedyPlusPlus(G, tau);
    std::vector<int> greedy_result =
        (greedy_full.size() > greedy_scope.size()) ? std::move(greedy_full)
                                                   : std::move(greedy_scope);

    // ── Step 3: connectivity-aware batch peeling on scope ────
    std::vector<int> best_result = greedy_result;
    int best_size = (int)greedy_result.size();

    // Reusable scratch.
    std::vector<bool> active(G.n, false);
    std::vector<int>  loc_deg(G.n, 0);

    // Work queue: each entry is the current active node-list for a state.
    std::queue<std::vector<int>> work_q;
    {
        std::vector<int> init_nodes;
        init_nodes.reserve(std::count(scope.begin(), scope.end(), true));
        for (int v = 0; v < G.n; ++v) if (scope[v]) init_nodes.push_back(v);
        if (!init_nodes.empty()) work_q.push(std::move(init_nodes));
    }

    while (!work_q.empty()) {
        std::vector<int> nodes = std::move(work_q.front());
        work_q.pop();
        int active_cnt = (int)nodes.size();
        if (active_cnt <= best_size) continue;

        // Activate.
        for (int v : nodes) active[v] = true;
        for (int v : nodes) {
            int d = 0;
            for (int u : G.adj[v]) if (active[u]) ++d;
            loc_deg[v] = d;
        }

        auto tryUpdate = [&]() {
            if (active_cnt <= best_size) return;
            int thr = floorPow(active_cnt, tau);
            for (int v : nodes)
                if (active[v] && loc_deg[v] < thr) return;
            if (!G.isConnected(active)) return;
            best_size = active_cnt;
            best_result.clear();
            for (int v : nodes) if (active[v]) best_result.push_back(v);
        };

        tryUpdate();

        // Batch-peel loop.  Each round: Tarjan AP once, remove every
        // non-AP below threshold, then detect splits and enqueue all
        // large-enough child CCs.
        bool changed = true;
        while (changed && active_cnt > best_size) {
            changed = false;
            int thr = floorPow(active_cnt, tau);

            int min_deg = INT_MAX;
            for (int v : nodes)
                if (active[v]) min_deg = std::min(min_deg, loc_deg[v]);
            if (min_deg >= thr) { tryUpdate(); break; }

            std::vector<bool> ap = G.articulationPoints(active);

            std::vector<int> to_remove;
            for (int v : nodes)
                if (active[v] && !ap[v] && loc_deg[v] < thr)
                    to_remove.push_back(v);

            if (to_remove.empty()) break;

            for (int v : to_remove) {
                active[v] = false;
                --active_cnt;
                for (int w : G.adj[v]) if (active[w]) --loc_deg[w];
                loc_deg[v] = 0;
            }
            changed = true;

            // Detect splits among the originally-listed nodes.
            std::vector<bool> vis(G.n, false);
            int largest_sz = 0;
            std::vector<int> largest_comp;

            for (int s0 : nodes) {
                if (!active[s0] || vis[s0]) continue;
                std::vector<int> comp;
                std::queue<int> bq;
                bq.push(s0); vis[s0] = true;
                while (!bq.empty()) {
                    int u = bq.front(); bq.pop();
                    comp.push_back(u);
                    for (int w : G.adj[u])
                        if (active[w] && !vis[w]) { vis[w] = true; bq.push(w); }
                }

                if ((int)comp.size() > largest_sz) {
                    // Previous "largest" becomes a sibling → enqueue if useful.
                    if (!largest_comp.empty() &&
                        (int)largest_comp.size() > best_size) {
                        work_q.push(std::move(largest_comp));
                    }
                    largest_sz   = (int)comp.size();
                    largest_comp = std::move(comp);
                } else if ((int)comp.size() > best_size) {
                    work_q.push(std::move(comp));
                }
            }

            if (largest_sz == 0) break;  // nothing left
            if (largest_sz < active_cnt) {
                // Split occurred: keep peeling the largest CC in-place.
                std::vector<bool> in_largest(G.n, false);
                for (int v : largest_comp) in_largest[v] = true;
                for (int v : nodes) {
                    if (active[v] && !in_largest[v]) active[v] = false;
                }
                active_cnt = largest_sz;
                // Recompute loc_deg within the largest CC only.
                for (int v : largest_comp) {
                    int d = 0;
                    for (int u : G.adj[v]) if (active[u]) ++d;
                    loc_deg[v] = d;
                }
                nodes = std::move(largest_comp);
            }
        }

        // Deactivate everything this state touched.
        for (int v : nodes) { active[v] = false; loc_deg[v] = 0; }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    AlgoResult res;
    res.nodes    = best_result;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    return res;
}
