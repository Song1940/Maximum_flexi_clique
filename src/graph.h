#pragma once
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <functional>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <vector>
#include <string>

struct Graph {
    int n = 0, m = 0;
    std::vector<std::vector<int>> adj;  // sorted adjacency lists

    Graph() = default;
    explicit Graph(int n) : n(n), adj(n) {}

    void addEdge(int u, int v) {
        assert(u >= 0 && u < n && v >= 0 && v < n && u != v);
        adj[u].push_back(v);
        adj[v].push_back(u);
        ++m;
    }

    // Remove duplicate edges and self-loops; sort adjacency lists.
    void finalize();

    bool hasEdge(int u, int v) const {
        return std::binary_search(adj[u].begin(), adj[u].end(), v);
    }

    int deg(int v) const { return (int)adj[v].size(); }

    // Bucket-based k-core decomposition, O(n + m).
    std::vector<int> coreNumbers() const;

    // Tarjan's articulation-point detection on the active subgraph.
    std::vector<bool> articulationPoints(const std::vector<bool>& active) const;

    // Returns true if removing v disconnects its active neighbors.
    bool isArticulationPoint(int v, const std::vector<bool>& active) const;

    // BFS distances from src within the active subgraph; unreachable = INT_MAX.
    std::vector<int> bfs(int src, const std::vector<bool>& active) const;

    // Returns true if all active nodes form a single connected component.
    bool isConnected(const std::vector<bool>& active) const;

    // Returns a bool mask of the largest connected component among active nodes.
    std::vector<bool> largestCC(const std::vector<bool>& active) const;

    // Load an edge-list file (whitespace-separated); re-indexes nodes to [0, n).
    static Graph loadFromFile(const std::string& path);
};

// floor(k ^ tau)
inline int floorPow(long long k, double tau) {
    if (k <= 0) return 0;
    return static_cast<int>(std::floor(std::pow(static_cast<double>(k), tau)));
}

// Minimum degree threshold needed to beat best_size.
// A flexi-clique that improves on best_size has size >= best_size+1, so every
// member needs internal degree >= floor((best_size+1)^tau); a node below that can
// belong to no improving flexi-clique and is removed (Rule 1 and the global degree
// prune). Hence theta = floor((best_size+1)^tau).
//
// This floor must be computed exactly: too small weakens pruning, too large
// over-prunes and could drop optimal solutions. We take the long-double floor
// and correct boundary rounding with the exact inverse test m <= x^tau iff
// m^(1/tau) <= x, so theta is accepted only when verifiably <= x^tau.
inline int computeTheta(int best_size, double tau) {
    const long double x = static_cast<long double>(best_size) + 1.0L;
    long long m = static_cast<long long>(std::floor(std::pow(x, static_cast<long double>(tau))));
    if (m < 0) m = 0;
    auto fits = [&](long long t) -> bool {            // t <= x^tau  <=>  t^(1/tau) <= x
        return t <= 0 || std::pow(static_cast<long double>(t),
                                  1.0L / static_cast<long double>(tau)) <= x;
    };
    while (fits(m + 1)) ++m;
    while (m > 0 && !fits(m)) --m;
    return static_cast<int>(m);
}

// Lemma 2: minimum number of nodes in a connected graph
// with minimum degree k and diameter L.
inline int nFunc(int k, int L) {
    if (L <= 0) return 1;
    if (k <= 0) return L + 1;
    if (L <= 2 || k == 1) return k + L;
    return k + L + 1 + (L / 3) * (k - 2);
}

// Returns true if nodes form a valid tau-flexi-clique in G.
bool isFlexiClique(const Graph& G, const std::vector<int>& nodes, double tau);
