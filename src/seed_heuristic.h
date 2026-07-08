#pragma once
#include "graph.h"
#include "npa.h"

// Multi-seed greedy heuristic for the EBA initial lower bound.
AlgoResult multiSeedGreedyPlusPlus(const Graph& G, double tau);
