#pragma once
#include "graph.h"
#include <vector>
#include <string>

struct AlgoResult {
    std::vector<int> nodes;
    double time_sec = 0.0;
    long long branches     = 0;  // only meaningful for EBA
    long long depth_prunes = 0;  // depth rule fires (eba_pivot only)
    bool timed_out = false;      // seeded anytime mode hit its deadline
};

// NPA: Naive Peeling Algorithm
// K-core seeding followed by articulation-point-preserving peeling.
AlgoResult runNPA(const Graph& G, double tau);
