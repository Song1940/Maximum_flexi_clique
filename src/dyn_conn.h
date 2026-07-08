#pragma once
// ---- Link-Cut Tree (1-indexed): link/cut/connected, O(log n) amortised ----
class LCT {
    struct Node {
        int ch[2], fa;
        bool rev;
    };
    std::vector<Node> nd;
    std::vector<int>  _buf;  // reusable stack for pushAll

    bool isRoot(int x) const {
        int p = nd[x].fa;
        return nd[p].ch[0] != x && nd[p].ch[1] != x;
    }

    void flip(int x) {
        std::swap(nd[x].ch[0], nd[x].ch[1]);
        nd[x].rev ^= 1;
    }

    void push(int x) {
        if (nd[x].rev) {
            if (nd[x].ch[0]) flip(nd[x].ch[0]);
            if (nd[x].ch[1]) flip(nd[x].ch[1]);
            nd[x].rev = 0;
        }
    }

    // Push lazy from splay-root down to x (iterative).
    void pushAll(int x) {
        _buf.clear();
        for (int y = x; !isRoot(y); y = nd[y].fa) _buf.push_back(y);
        // nd[_buf.back()].fa is the splay-root (if _buf non-empty), else x itself.
        int sr = _buf.empty() ? x : nd[_buf.back()].fa;
        push(sr);
        for (int i = (int)_buf.size() - 1; i >= 0; --i) push(_buf[i]);
    }

    void rot(int x) {
        int y = nd[x].fa, z = nd[y].fa;
        int k = (nd[y].ch[1] == x);
        if (!isRoot(y)) nd[z].ch[nd[z].ch[1] == y] = x;
        nd[x].fa = z;
        nd[y].ch[k] = nd[x].ch[k ^ 1];
        if (nd[x].ch[k ^ 1]) nd[nd[x].ch[k ^ 1]].fa = y;
        nd[x].ch[k ^ 1] = y;
        nd[y].fa = x;
    }

    void splay(int x) {
        pushAll(x);
        while (!isRoot(x)) {
            int y = nd[x].fa;
            if (!isRoot(y)) {
                int z = nd[y].fa;
                rot((nd[z].ch[0] == y) == (nd[y].ch[0] == x) ? y : x);
            }
            rot(x);
        }
    }

    void access(int x) {
        for (int y = x, z = 0; y; y = nd[y].fa) {
            splay(y);
            nd[y].ch[1] = z;
            z = y;
        }
        splay(x);
    }

    void makeRoot(int x) {
        access(x);
        flip(x);
    }

    int findRoot(int x) {
        access(x);
        while (nd[x].ch[0]) { push(x); x = nd[x].ch[0]; }
        splay(x);  // amortize
        return x;
    }

public:
    explicit LCT(int n) : nd(n + 1, {{0, 0}, 0, false}), _buf() {
        _buf.reserve(64);
    }

    // Are u and v in the same tree?
    bool connected(int u, int v) {
        makeRoot(u);
        return findRoot(v) == u;
    }

    // Add edge u–v.  Caller must ensure u and v are NOT already connected.
    void link(int u, int v) {
        makeRoot(u);
        nd[u].fa = v;
    }

    // Remove edge u–v.  Caller must ensure the edge exists.
    void cut(int u, int v) {
        makeRoot(u);
        access(v);
        // u is now the leftmost node in v's splay tree (direct left child).
        nd[v].ch[0] = 0;
        nd[u].fa    = 0;
    }
};
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>

// Dynamic connectivity for the FPA peeling.
//
// A spanning forest of the current graph is maintained in a Link-Cut Tree (the
// dynamic-connectivity structure cited by the paper). Connectivity queries and
// the component count are answered exactly. When a tree edge is deleted, a
// replacement non-tree edge crossing the cut is sought by enumerating the
// smaller of the two spanning-tree sides (dual BFS over the tree adjacency) and
// scanning its incident non-tree edges. Since the component count is a property
// of the graph, the result is independent of which replacement is chosen.
class FullyDynamicConnectivity {
    int n;
    long long comp_count;
    LCT forest;                                              // node v -> v+1
    std::unordered_map<int, std::unordered_set<int>> tadj;   // spanning-tree adjacency
    std::unordered_map<int, std::unordered_set<int>> nadj;   // non-tree adjacency
    std::unordered_set<long long> present;

    static long long key(int u, int v) { if (u > v) std::swap(u, v); return ((long long)u << 32) | (unsigned)v; }
    bool fc(int u, int v) { return forest.connected(u + 1, v + 1); }

    void tlink(int u, int v) { forest.link(u + 1, v + 1); tadj[u].insert(v); tadj[v].insert(u); }
    void tcut (int u, int v) { forest.cut (u + 1, v + 1); tadj[u].erase(v);  tadj[v].erase(u);  }

    // Enumerate the smaller spanning-tree side of a just-cut edge (u,v) by
    // expanding both sides one node at a time; the side that finishes first is
    // the smaller one.
    void enumerateSmaller(int u, int v, std::vector<int>& out) {
        std::unordered_set<int> vu{u}, vv{v};
        std::queue<int> qu, qv; qu.push(u); qv.push(v);
        while (!qu.empty() && !qv.empty()) {
            int x = qu.front(); qu.pop();
            auto it = tadj.find(x);
            if (it != tadj.end()) for (int w : it->second) if (vu.insert(w).second) qu.push(w);
            int y = qv.front(); qv.pop();
            auto jt = tadj.find(y);
            if (jt != tadj.end()) for (int w : jt->second) if (vv.insert(w).second) qv.push(w);
        }
        if (qu.empty()) out.assign(vu.begin(), vu.end());
        else            out.assign(vv.begin(), vv.end());
    }

public:
    explicit FullyDynamicConnectivity(int n_) : n(n_), comp_count(n_), forest(n_) {}

    bool connected(int u, int v) { return fc(u, v); }
    long long get_component_count() const { return comp_count; }

    void insert_edge(int u, int v) {
        if (u == v) return;
        long long k = key(u, v);
        if (!present.insert(k).second) return;
        if (!fc(u, v)) { tlink(u, v); --comp_count; }
        else           { nadj[u].insert(v); nadj[v].insert(u); }
    }

    void delete_edge(int u, int v) {
        if (u == v) return;
        long long k = key(u, v);
        if (!present.erase(k)) return;
        auto itu = nadj.find(u);
        if (itu != nadj.end() && itu->second.erase(v)) { nadj[v].erase(u); return; }  // non-tree
        tcut(u, v);
        std::vector<int> side;
        enumerateSmaller(u, v, side);
        std::unordered_set<int> inSide(side.begin(), side.end());
        int rx = -1, ry = -1;
        for (int x : side) {
            auto it = nadj.find(x);
            if (it == nadj.end()) continue;
            for (int y : it->second) if (!inSide.count(y)) { rx = x; ry = y; break; }
            if (rx != -1) break;
        }
        if (rx != -1) { nadj[rx].erase(ry); nadj[ry].erase(rx); tlink(rx, ry); }
        else          { ++comp_count; }
    }
};
