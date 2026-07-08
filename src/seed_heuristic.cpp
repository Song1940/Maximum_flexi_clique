#include "seed_heuristic.h"
#include "fpa.h"
#include "mgpp.h"
#include "graph.h"
#include <algorithm>
#include <chrono>
#include <numeric>
#include <vector>

AlgoResult multiSeedGreedyPlusPlus(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    AlgoResult base = runFPA(G, tau);
    std::vector<int> best = base.nodes;
    int best_size = (int)best.size();

    if (G.n == 0) {
        auto t1 = std::chrono::high_resolution_clock::now();
        AlgoResult r; r.nodes = best;
        r.time_sec = std::chrono::duration<double>(t1 - t0).count();
        return r;
    }

    // Vertices ordered by degree, descending.
    std::vector<int> order(G.n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b){ return G.adj[a].size() > G.adj[b].size(); });

    const int MAX_SEEDS = 40;
    const int BALL_CAP  = 5000;
    const int RADIUS    = 2;

    std::vector<bool> ball_mask(G.n, false);
    std::vector<bool> used_as_seed(G.n, false);
    std::vector<int>  ball_nodes;
    ball_nodes.reserve(BALL_CAP);

    int seeds_run = 0;
    for (int si = 0; si < G.n && seeds_run < MAX_SEEDS; ++si) {
        int s = order[si];
        if (used_as_seed[s]) continue;
        if ((int)G.adj[s].size() < 2) break;   // degree too low to matter

        // BFS 2-hop ball, capped.
        ball_nodes.clear();
        ball_mask[s] = true;
        ball_nodes.push_back(s);

        int layer_start = 0;
        for (int r = 0; r < RADIUS; ++r) {
            int layer_end = (int)ball_nodes.size();
            bool capped = false;
            for (int i = layer_start; i < layer_end && !capped; ++i) {
                int u = ball_nodes[i];
                for (int w : G.adj[u]) {
                    if (!ball_mask[w]) {
                        ball_mask[w] = true;
                        ball_nodes.push_back(w);
                        if ((int)ball_nodes.size() >= BALL_CAP) { capped = true; break; }
                    }
                }
            }
            layer_start = layer_end;
            if (capped) break;
        }

        // Mark the ball as used so neighbouring vertices are skipped as later seeds.
        for (int u : ball_nodes) used_as_seed[u] = true;

        if ((int)ball_nodes.size() > best_size) {
            std::vector<int> local = modifiedGreedyPlusPlus(G, tau, &ball_mask);
            if ((int)local.size() > best_size) {
                best      = std::move(local);
                best_size = (int)best.size();
            }
        }

        // Clear the mask for the next iteration.
        for (int u : ball_nodes) ball_mask[u] = false;
        ++seeds_run;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    res.nodes    = best;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    return res;
}
