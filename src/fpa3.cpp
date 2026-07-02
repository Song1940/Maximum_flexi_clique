#include "fpa3.h"
#include "fpa.h"
#include "graph.h"
#include <algorithm>
#include <chrono>
#include <numeric>
#include <vector>

// ============================================================
// fpa3.cpp — FPA with community-aware multi-seed MGPP.
//
// Motivation: on LFR-style synthetic graphs with near-uniform
// degree, the k-core seed used by the original FPA is a poor
// proxy for the true densest community. MGPP on the full graph
// also struggles because every region looks similar.  But the
// true maximum flexi-clique lives inside a single community,
// whose diameter is usually ≤2. So expanding a 2-hop ball
// around a high-degree vertex captures that community plus
// some fringe, and MGPP restricted to that ball converges on it.
//
// Algorithm
// ---------
// 1) Run original runFPA as a baseline.
// 2) Order vertices by degree, descending.
// 3) For at most MAX_SEEDS unused top-degree vertices s:
//      - Expand a 2-hop ball around s, capped at BALL_CAP.
//      - Mark every node in the ball as "used" so later seeds
//        pick a different region.
//      - If the ball is larger than the current best, run MGPP
//        restricted to the ball; accept any improvement.
// ============================================================

AlgoResult runFPA3(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Step 1: baseline FPA.
    AlgoResult base = runFPA(G, tau);
    std::vector<int> best = base.nodes;
    int best_size = (int)best.size();

    if (G.n == 0) {
        auto t1 = std::chrono::high_resolution_clock::now();
        AlgoResult r; r.nodes = best;
        r.time_sec = std::chrono::duration<double>(t1 - t0).count();
        return r;
    }

    // Step 2: vertices ordered by degree descending.
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

        // Mark ball as used so neighbouring vertices are skipped as later seeds.
        for (int u : ball_nodes) used_as_seed[u] = true;

        if ((int)ball_nodes.size() > best_size) {
            std::vector<int> local = modifiedGreedyPlusPlus(G, tau, &ball_mask);
            if ((int)local.size() > best_size) {
                best      = std::move(local);
                best_size = (int)best.size();
            }
        }

        // Clear mask for next iteration.
        for (int u : ball_nodes) ball_mask[u] = false;
        ++seeds_run;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    res.nodes    = best;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    return res;
}
