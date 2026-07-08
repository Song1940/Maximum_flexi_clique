#pragma once
#include "graph.h"
#include <vector>

// Greedy peeling heuristic for the EBA initial lower bound; if active_mask is
// non-null, only masked nodes are considered.
std::vector<int> modifiedGreedyPlusPlus(const Graph& G, double tau,
                                        const std::vector<bool>* active_mask = nullptr);
