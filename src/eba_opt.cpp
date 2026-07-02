#include "eba.h"
#include "fpa.h"
#include "fpa3.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <queue>
#include <unordered_set>


struct EBAState {
    const Graph& G;
    double tau;

    std::vector<int> best;
    int best_size = 0;
    int theta = 0;

    std::vector<bool> g_active;
    std::vector<int>  g_adj_deg;
    std::vector<int> local_adj_deg;
    std::vector<int> dist_S;

    long long branches = 0;
    int force_seed = -1;

    std::chrono::steady_clock::time_point deadline;
    bool has_deadline = false;
    bool timed_out = false;

    std::vector<bool> bfs_scope_arr;
    std::vector<int>  bfs_dist_arr;
    std::vector<int>  bfs_queue;
    std::vector<int>  bfs_scope_touched;
    std::vector<int>  bfs_dist_touched;

    mutable std::vector<bool> removed_sim_arr;
    mutable std::vector<int>  sim_deg_arr;   // INT_MIN = uninitialized
    mutable std::vector<int>  sim_d_touched;
    mutable std::vector<int>  sim_r_touched;
    mutable std::vector<int>  sim_bfs;

    std::vector<bool> cun_arr;

    // OPT-1: incremental S-membership and S-internal degree.
    std::vector<bool> in_S;
    std::vector<int>  deg_in_S;

    EBAState(const Graph& G, double tau)
        : G(G), tau(tau),
          g_active(G.n, true), g_adj_deg(G.n, 0),
          local_adj_deg(G.n, 0), dist_S(G.n, INT_MAX),
          bfs_scope_arr(G.n, false), bfs_dist_arr(G.n, INT_MAX),
          removed_sim_arr(G.n, false), sim_deg_arr(G.n, INT_MIN),
          cun_arr(G.n, false),
          in_S(G.n, false), deg_in_S(G.n, 0) {
        for (int v = 0; v < G.n; ++v)
            g_adj_deg[v] = (int)G.adj[v].size();
    }

    void updateTheta() { theta = computeTheta(best_size, tau); }

    std::vector<int> globalDegreePrune() {
        std::vector<int> removed;
#ifndef NO_RULE6
        std::queue<int> q;
        for (int v = 0; v < G.n; ++v)
            if (g_active[v] && g_adj_deg[v] < theta) q.push(v);
        while (!q.empty()) {
            int v = q.front(); q.pop();
            if (!g_active[v]) continue;
            g_active[v] = false;
            removed.push_back(v);
            for (int u : G.adj[v]) {
                if (g_active[u] && --g_adj_deg[u] < theta) q.push(u);
            }
        }
#endif
        return removed;
    }

    std::vector<int> computeFollowers(int v,
                                      const std::vector<bool>& in_scope) const {
        std::vector<int> followers;
#ifndef NO_RULE5
        sim_d_touched.clear();
        sim_r_touched.clear();
        sim_bfs.clear();

        auto initDeg = [&](int u) -> int& {
            if (sim_deg_arr[u] == INT_MIN) {
                sim_deg_arr[u] = local_adj_deg[u];
                sim_d_touched.push_back(u);
            }
            return sim_deg_arr[u];
        };

        auto markRm = [&](int u) {
            if (!removed_sim_arr[u]) {
                removed_sim_arr[u] = true;
                sim_r_touched.push_back(u);
            }
        };

        markRm(v);
        for (int u : G.adj[v]) {
            if (!in_scope[u] || removed_sim_arr[u]) continue;
            if (--initDeg(u) < theta) { markRm(u); sim_bfs.push_back(u); }
        }

        for (int qi = 0; qi < (int)sim_bfs.size(); ++qi) {
            int u = sim_bfs[qi];
            followers.push_back(u);
            for (int w : G.adj[u]) {
                if (!in_scope[w] || removed_sim_arr[w]) continue;
                if (--initDeg(w) < theta) { markRm(w); sim_bfs.push_back(w); }
            }
        }

        for (int u : sim_d_touched) sim_deg_arr[u]     = INT_MIN;
        for (int u : sim_r_touched) removed_sim_arr[u] = false;
#else
        (void)v; (void)in_scope;
#endif
        return followers;
    }

    void tryUpdate(const std::vector<int>& S) {
        if ((int)S.size() <= best_size) return;
        int k = (int)S.size();
        if (k >= 2) {
            int threshold = floorPow(k, tau);
            for (int v : S)
                if (deg_in_S[v] < threshold) return;
        }
        best      = S;
        best_size = k;
        updateTheta();
        globalDegreePrune();
    }

    void recurse(std::vector<int>&  S,
                 std::vector<int>   Cr,
                 std::vector<int>   Cun,
                 std::vector<bool>& in_D,
                 std::vector<bool>& in_scope) {
        ++branches;

        if (has_deadline && (branches & 4095) == 0 &&
            std::chrono::steady_clock::now() >= deadline) timed_out = true;
        if (timed_out) return;
        int min_adj_S = INT_MAX;
#if !defined(NO_RULE1) || !defined(NO_RULE2)
        if (!S.empty())
            for (int v : S) min_adj_S = std::min(min_adj_S, local_adj_deg[v]);
#endif

#ifndef NO_RULE1
        if (!S.empty() && min_adj_S < theta) return;
#endif

#ifndef NO_RULE2
        if (!S.empty()) {
            int max_possible = (int)std::floor(std::pow(min_adj_S + 1.0, 1.0 / tau));
            if ((int)S.size() >= max_possible) return;
        }
#endif

        if ((int)S.size() + (int)Cr.size() + (int)Cun.size() <= best_size) return;
        std::vector<int> forced_root;
        if (S.empty() && force_seed >= 0) forced_root.push_back(force_seed);
        std::vector<int>& R =
            (S.empty() && force_seed >= 0) ? forced_root : (S.empty() ? Cun : Cr);
#ifndef NO_DEGREE_ORDER
        std::sort(R.begin(), R.end(),
                  [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
#endif
        for (int u : Cun) cun_arr[u] = true;
        std::vector<std::pair<int,int>> branch_deg_undo;
        std::vector<std::pair<int,int>> branch_dist_undo;
        std::vector<int> newly_excluded, dist_pruned;
        std::vector<int> new_Cr, new_Cun, v_cn;

        std::vector<std::pair<int,int>> sibling_deg_undo;
        std::vector<int> prev_followers;
        int prev_sibling = -1;

        for (int idx = 0; idx < (int)R.size(); ++idx) {
            int v = R[idx];
            if (in_D[v] || !g_active[v]) continue;

            branch_deg_undo.clear();
            branch_dist_undo.clear();
            newly_excluded.clear();
            dist_pruned.clear();
            new_Cr.clear(); new_Cun.clear(); v_cn.clear();

            // Exclude followers of the previous sibling before branching on v.
            for (int f : prev_followers) {
                if (in_D[f] || !in_scope[f]) continue;
                in_D[f] = true; in_scope[f] = false;
                newly_excluded.push_back(f);
                for (int u : G.adj[f]) {
                    if (in_scope[u]) {
                        branch_deg_undo.push_back({u, local_adj_deg[u]});
                        --local_adj_deg[u];
                    }
                }
            }

            int new_scope = (int)S.size() + (int)Cr.size() + (int)Cun.size();
            if (new_scope <= best_size) {
                for (int e : newly_excluded) { in_scope[e] = true; in_D[e] = false; }
                for (auto& [nd, old] : branch_deg_undo) local_adj_deg[nd] = old;
                break;
            }

            S.push_back(v);
            in_scope[v] = true;

            // OPT-1: maintain deg_in_S incrementally.
            in_S[v] = true;
            for (int u : G.adj[v])
                if (in_S[u]) { ++deg_in_S[v]; ++deg_in_S[u]; }

            // Build v's Cun-neighbors; temporarily unmark in cun_arr so they
            // are excluded from new_Cun and promoted to new_Cr.
            // FIX: also exclude nodes already in S (!in_S[u]) to prevent stale
            // cun_arr values (left over from parent-level restores) from
            // re-introducing S members into Cr, which would cause S duplicates.
            for (int u : G.adj[v]) {
                if (cun_arr[u] && !in_D[u] && g_active[u] && !in_S[u]) {
                    v_cn.push_back(u);
                    cun_arr[u] = false;
                }
            }

            for (int u : Cr)  { if (u != v && !in_D[u] && g_active[u] && !in_S[u]) new_Cr.push_back(u); }
            for (int u : v_cn) new_Cr.push_back(u);
            for (int u : Cun) { if (u != v && cun_arr[u] && !in_D[u] && g_active[u] && !in_S[u]) new_Cun.push_back(u); }

            for (int u : v_cn) cun_arr[u] = true;  // restore before child call

#if !defined(NO_RULE3) || !defined(NO_RULE4)
            static constexpr int BFS_SCOPE_LIMIT = 30000;
            if ((int)S.size() >= 2 ||
                (int)(new_Cr.size() + new_Cun.size()) <= BFS_SCOPE_LIMIT) {

                auto markScope = [&](int u) {
                    if (!bfs_scope_arr[u]) {
                        bfs_scope_arr[u] = true;
                        bfs_scope_touched.push_back(u);
                    }
                };
                for (int u : S)       markScope(u);
                for (int u : new_Cr)  markScope(u);
                for (int u : new_Cun) markScope(u);

                bfs_queue.clear();
                bfs_dist_arr[v] = 0;
                bfs_dist_touched.push_back(v);
                bfs_queue.push_back(v);

                for (int qi = 0; qi < (int)bfs_queue.size(); ++qi) {
                    int u = bfs_queue[qi];
                    for (int w : G.adj[u]) {
                        if (bfs_scope_arr[w] && bfs_dist_arr[w] == INT_MAX) {
                            bfs_dist_arr[w] = bfs_dist_arr[u] + 1;
                            bfs_dist_touched.push_back(w);
                            bfs_queue.push_back(w);
                        }
                    }
                }

                for (int u : new_Cr) {
                    int d = bfs_dist_arr[u];
                    if (d != INT_MAX && d < dist_S[u]) {
                        branch_dist_undo.push_back({u, dist_S[u]});
                        dist_S[u] = d;
                    }
                }
                for (int u : new_Cun) {
                    int d = bfs_dist_arr[u];
                    if (d != INT_MAX && d < dist_S[u]) {
                        branch_dist_undo.push_back({u, dist_S[u]});
                        dist_S[u] = d;
                    }
                }

                for (int u : bfs_dist_touched)  bfs_dist_arr[u]  = INT_MAX;
                for (int u : bfs_scope_touched) bfs_scope_arr[u] = false;
                bfs_dist_touched.clear();
                bfs_scope_touched.clear();
            }
#endif

            [[maybe_unused]] int H_size =
                (int)S.size() + (int)new_Cr.size() + (int)new_Cun.size();

            [[maybe_unused]] int rule4_min_adj = INT_MAX;
#ifndef NO_RULE4
            if (!S.empty())
                for (int s : S) rule4_min_adj = std::min(rule4_min_adj, local_adj_deg[s]);
#endif

            auto filterByDist = [&](std::vector<int>& cands) {
                cands.erase(std::remove_if(cands.begin(), cands.end(), [&](int u) {
                    int du = dist_S[u];
                    if (du == INT_MAX) return false;
#ifndef NO_RULE3
                    if (nFunc(theta, du) > H_size) { dist_pruned.push_back(u); return true; }
#endif
#ifndef NO_RULE4
                    if (rule4_min_adj < floorPow(nFunc(theta, du), tau)) {
                        dist_pruned.push_back(u); return true;
                    }
#endif
                    return false;
                }), cands.end());
            };

#if !defined(NO_RULE3) || !defined(NO_RULE4)
            filterByDist(new_Cr);
            filterByDist(new_Cun);
            for (int u : dist_pruned) {
                if (!in_D[u]) {
                    in_D[u] = true; in_scope[u] = false;
                    for (int w : G.adj[u]) {
                        if (in_scope[w]) {
                            branch_deg_undo.push_back({w, local_adj_deg[w]});
                            --local_adj_deg[w];
                        }
                    }
                }
            }
#endif

            tryUpdate(S);

#ifndef NO_DEGREE_ORDER
            std::sort(new_Cr.begin(), new_Cr.end(),
                      [&](int a, int b){ return local_adj_deg[a] < local_adj_deg[b]; });
#endif

            // Rule 1 pre-check: avoid entering child branch that will prune immediately.
#ifndef NO_RULE1
            if (!S.empty()) {
                int child_min = INT_MAX;
                for (int u : S) child_min = std::min(child_min, local_adj_deg[u]);
                if (child_min < theta) goto BACKTRACK;
            }
#endif

            if ((int)S.size() + (int)new_Cr.size() + (int)new_Cun.size() > best_size) {
                recurse(S, new_Cr, new_Cun, in_D, in_scope);
                // Restore cun_arr for child's Cun after child has reset it.
                for (int u : new_Cun) cun_arr[u] = true;
            }

            [[maybe_unused]] BACKTRACK:
            S.pop_back();
            in_scope[v] = false;

            // OPT-1: undo deg_in_S.
            in_S[v] = false;
            for (int u : G.adj[v])
                if (in_S[u]) { --deg_in_S[v]; --deg_in_S[u]; }

            for (auto& [nd, old] : branch_dist_undo) dist_S[nd]        = old;
            for (auto& [nd, old] : branch_deg_undo)  local_adj_deg[nd] = old;
            for (int u : dist_pruned)    if (in_D[u]) { in_D[u] = false; in_scope[u] = true; }
            for (int e : newly_excluded) if (in_D[e]) { in_D[e] = false; in_scope[e] = true; }

#ifndef NO_RULE5
            prev_followers = computeFollowers(v, in_scope);
#else
            prev_followers.clear();
#endif
            prev_sibling = v;

            in_D[v] = true; in_scope[v] = false;
            for (int u : G.adj[v]) {
                if (in_scope[u]) {
                    sibling_deg_undo.push_back({u, local_adj_deg[u]});
                    --local_adj_deg[u];
                }
            }
        }

        for (auto it = sibling_deg_undo.rbegin(); it != sibling_deg_undo.rend(); ++it)
            local_adj_deg[it->first] = it->second;

        if (prev_sibling != -1) {
            for (int u : R)             if (in_D[u]) { in_D[u] = false; in_scope[u] = true; }
            for (int f : prev_followers) if (in_D[f]) { in_D[f] = false; in_scope[f] = true; }
        }

        for (int u : Cun) cun_arr[u] = false;
    }
};

AlgoResult runEBA(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();
    
    AlgoResult fpa_res = runFPA3(G, tau);
    AlgoResult npa_res = runNPA(G, tau);

    std::vector<int> init_sol = fpa_res.nodes;
    if (npa_res.nodes.size() > init_sol.size()) init_sol = npa_res.nodes;

    EBAState state(G, tau);
    state.best      = init_sol;
    state.best_size = (int)init_sol.size();
    state.updateTheta();

    state.globalDegreePrune();

    {
        std::vector<int> mgpp_pruned =
            modifiedGreedyPlusPlus(G, tau, &state.g_active);
        if ((int)mgpp_pruned.size() > state.best_size) {
            state.best      = std::move(mgpp_pruned);
            state.best_size = (int)state.best.size();
            state.updateTheta();
            state.globalDegreePrune();
        }
    }

#ifdef NO_RULE6
    state.best_size = 0;
    state.best.clear();
    state.updateTheta();
#endif

    std::vector<int> V_prime;
    for (int v = 0; v < G.n; ++v)
        if (state.g_active[v]) V_prime.push_back(v);

    if (V_prime.empty()) {
        AlgoResult res;
        res.nodes    = init_sol;
        res.time_sec = fpa_res.time_sec;
        return res;
    }

    for (int v : V_prime) {
        int d = 0;
        for (int u : G.adj[v])
            if (state.g_active[u]) ++d;
        state.local_adj_deg[v] = d;
    }

    std::vector<int>  S, Cr, Cun = V_prime;
    std::vector<bool> in_D(G.n, false), in_scope(G.n, false);
    for (int v : V_prime) in_scope[v] = true;
    std::fill(state.dist_S.begin(), state.dist_S.end(), INT_MAX);

    state.recurse(S, Cr, Cun, in_D, in_scope);

    auto t1 = std::chrono::high_resolution_clock::now();

    AlgoResult res;
    // Pick the better of EBA result and heuristic initial solution.
    // (Relevant for NO_RULE6 where search starts with best_size=0.)
    if (state.best.size() >= init_sol.size())
        res.nodes = state.best.empty() ? init_sol : state.best;
    else
        res.nodes = init_sol;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    res.branches = state.branches;
    return res;
}

AlgoResult runEBASeeded(const Graph& region, int seed, double tau,
                        double timeout_sec) {
    auto t0 = std::chrono::high_resolution_clock::now();
    AlgoResult res;
    if (seed < 0 || seed >= region.n) { res.nodes = {}; return res; }

    EBAState state(region, tau);
    state.best       = { seed };          // singleton is always valid
    state.best_size  = 1;
    state.force_seed = seed;
    state.updateTheta();
    if (timeout_sec > 0) {
        state.has_deadline = true;
        state.deadline = std::chrono::steady_clock::now() +
            std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(timeout_sec));
    }

    for (int v = 0; v < region.n; ++v) {
        int d = 0;
        for (int u : region.adj[v]) if (state.g_active[u]) ++d;
        state.local_adj_deg[v] = d;
    }

    std::vector<int>  S, Cr, Cun;
    std::vector<bool> in_D(region.n, false), in_scope(region.n, false);
    for (int v = 0; v < region.n; ++v) { Cun.push_back(v); in_scope[v] = true; }
    std::fill(state.dist_S.begin(), state.dist_S.end(), INT_MAX);

    state.recurse(S, Cr, Cun, in_D, in_scope);

    auto t1 = std::chrono::high_resolution_clock::now();
    res.nodes     = state.best.empty() ? std::vector<int>{seed} : state.best;
    res.time_sec  = std::chrono::duration<double>(t1 - t0).count();
    res.branches  = state.branches;
    res.timed_out = state.timed_out;
    return res;
}
