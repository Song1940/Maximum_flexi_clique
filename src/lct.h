#pragma once
#include <vector>

// Link-Cut Tree (1-indexed).
// Supports link(u,v), cut(u,v), connected(u,v) in O(log n) amortized.
// Based on access / makeRoot / findRoot idiom.
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
