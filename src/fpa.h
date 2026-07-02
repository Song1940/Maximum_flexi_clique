#pragma once
#include "graph.h"
#include "npa.h"

// FPA: Flexi-Prune Algorithm
// K-core seed + Modified Greedy++ initialization + connectivity-aware peeling.
AlgoResult runFPA(const Graph& G, double tau);

// Modified Greedy++ heuristic used internally by FPA and EBA.
// If active_mask is non-null, only nodes with active_mask[v]=true are considered.
std::vector<int> modifiedGreedyPlusPlus(const Graph& G, double tau,
                                        const std::vector<bool>* active_mask = nullptr);
