#include "graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

void Graph::finalize() {
    m = 0;
    for (int v = 0; v < n; ++v) {
        auto& a = adj[v];
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
        auto it = std::find(a.begin(), a.end(), v);
        if (it != a.end()) a.erase(it);
    }
    for (int v = 0; v < n; ++v) m += (int)adj[v].size();
    m /= 2;
}

// Bucket-based k-core decomposition (O(n + m)).
std::vector<int> Graph::coreNumbers() const {
    std::vector<int> core(n, 0);
    std::vector<int> deg(n);
    for (int v = 0; v < n; ++v) deg[v] = (int)adj[v].size();

    int max_d = *std::max_element(deg.begin(), deg.end());

    std::vector<int> bin(max_d + 1, 0);
    for (int v = 0; v < n; ++v) bin[deg[v]]++;

    int start = 0;
    std::vector<int> pos(max_d + 1);
    for (int d = 0; d <= max_d; ++d) { pos[d] = start; start += bin[d]; }

    std::vector<int> vert(n), pos_v(n);
    for (int v = 0; v < n; ++v) {
        pos_v[v]         = pos[deg[v]];
        vert[pos_v[v]]   = v;
        pos[deg[v]]++;
    }
    for (int d = max_d; d > 0; --d) pos[d] = pos[d - 1];
    pos[0] = 0;

    for (int i = 0; i < n; ++i) {
        int v = vert[i];
        core[v] = deg[v];
        for (int u : adj[v]) {
            if (core[u] > core[v]) {
                int du = deg[u], pu = pos_v[u], pw = pos[du], w = vert[pw];
                if (u != w) {
                    pos_v[u] = pw; vert[pu] = w;
                    pos_v[w] = pu; vert[pw] = u;
                }
                pos[du]++;
                deg[u]--;
            }
        }
    }
    return core;
}

// Tarjan's articulation-point algorithm on the active subgraph.
// Iterative (explicit stack) to avoid O(n) call-stack blow-up at 1M+ nodes.
std::vector<bool> Graph::articulationPoints(const std::vector<bool>& active) const {
    std::vector<bool> ap(n, false);
    std::vector<int>  disc(n, -1), low(n, 0), parent(n, -1);
    int timer = 0;

    // Stack entry: (u, next-neighbour-index into adj[u])
    std::vector<std::pair<int,int>> stk;
    stk.reserve(n);

    for (int start = 0; start < n; ++start) {
        if (!active[start] || disc[start] != -1) continue;

        disc[start] = low[start] = timer++;
        parent[start] = -1;
        int root_children = 0;
        stk.emplace_back(start, 0);

        while (!stk.empty()) {
            int u = stk.back().first;
            int& i = stk.back().second;  // ref is safe: reserve(n) => no realloc

            bool descended = false;
            while (i < (int)adj[u].size()) {
                int v = adj[u][i++];
                if (!active[v]) continue;
                if (disc[v] == -1) {
                    // Tree edge: descend
                    disc[v] = low[v] = timer++;
                    parent[v] = u;
                    if (u == start) ++root_children;
                    stk.emplace_back(v, 0);
                    descended = true;
                    break;
                } else if (v != parent[u]) {
                    // Back edge
                    low[u] = std::min(low[u], disc[v]);
                }
            }

            if (!descended) {
                // Post-order: propagate low to parent, mark AP if cut-vertex.
                int u_disc = disc[u], u_low = low[u];
                stk.pop_back();
                if (!stk.empty()) {
                    int p = stk.back().first;
                    low[p] = std::min(low[p], u_low);
                    if (parent[p] != -1 && u_low >= disc[p]) ap[p] = true;
                    (void)u_disc;
                }
            }
        }

        // Root is AP iff it has ≥2 DFS-tree children.
        if (root_children > 1) ap[start] = true;
    }
    return ap;
}

// Check if removing v disconnects its active neighbors.
bool Graph::isArticulationPoint(int v, const std::vector<bool>& active) const {
    std::vector<int> neighbors;
    for (int u : adj[v])
        if (active[u]) neighbors.push_back(u);
    if (neighbors.size() <= 1) return false;

    std::vector<bool> visited(n, false);
    visited[v] = true;
    std::queue<int> q;
    q.push(neighbors[0]);
    visited[neighbors[0]] = true;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int w : adj[u])
            if (active[w] && !visited[w]) { visited[w] = true; q.push(w); }
    }
    for (int u : neighbors)
        if (!visited[u]) return true;
    return false;
}

std::vector<int> Graph::bfs(int src, const std::vector<bool>& active) const {
    std::vector<int> dist(n, INT_MAX);
    if (!active[src]) return dist;
    dist[src] = 0;
    std::queue<int> q;
    q.push(src);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u])
            if (active[v] && dist[v] == INT_MAX) { dist[v] = dist[u] + 1; q.push(v); }
    }
    return dist;
}

bool Graph::isConnected(const std::vector<bool>& active) const {
    int start = -1, cnt = 0;
    for (int v = 0; v < n; ++v)
        if (active[v]) { ++cnt; if (start == -1) start = v; }
    if (cnt <= 1) return true;

    int visited = 0;
    std::vector<bool> seen(n, false);
    std::queue<int> q;
    q.push(start); seen[start] = true; ++visited;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int w : adj[u])
            if (active[w] && !seen[w]) { seen[w] = true; ++visited; q.push(w); }
    }
    return visited == cnt;
}

std::vector<bool> Graph::largestCC(const std::vector<bool>& active) const {
    std::vector<bool> in_lcc(n, false), visited(n, false);
    int best_size = 0;
    std::vector<int> best_comp;

    for (int s = 0; s < n; ++s) {
        if (!active[s] || visited[s]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(s); visited[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int w : adj[u])
                if (active[w] && !visited[w]) { visited[w] = true; q.push(w); }
        }
        if ((int)comp.size() > best_size) { best_size = (int)comp.size(); best_comp = comp; }
    }

    for (int v : best_comp) in_lcc[v] = true;
    return in_lcc;
}

Graph Graph::loadFromFile(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot open: " + path);

    std::unordered_map<int,int> id_map;
    std::vector<std::pair<int,int>> edges;
    int next_id = 0;

    auto getID = [&](int raw) -> int {
        auto [it, inserted] = id_map.emplace(raw, next_id);
        if (inserted) ++next_id;
        return it->second;
    };

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int u, v;
        if (!(ss >> u >> v) || u == v) continue;
        edges.push_back({getID(u), getID(v)});
    }

    Graph G(next_id);
    for (auto [u, v] : edges) G.addEdge(u, v);
    G.finalize();
    std::cout << "Loaded graph: n=" << G.n << " m=" << G.m << "\n";
    return G;
}

bool isFlexiClique(const Graph& G, const std::vector<int>& nodes, double tau) {
    if (nodes.empty()) return false;
    int k = (int)nodes.size();
    if (k == 1) return true;
    int threshold = floorPow(k, tau);

    std::unordered_set<int> nodeSet(nodes.begin(), nodes.end());
    for (int v : nodes) {
        int cnt = 0;
        for (int u : G.adj[v])
            if (nodeSet.count(u)) ++cnt;
        if (cnt < threshold) return false;
    }

    std::vector<bool> active(G.n, false);
    for (int v : nodes) active[v] = true;
    return G.isConnected(active);
}
