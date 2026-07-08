#include "mgpp.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <queue>
#include <utility>

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
