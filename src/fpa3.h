#pragma once
#include "graph.h"
#include "npa.h"

// FPA3: FPA + community-aware multi-seed MGPP probing.
// Designed for LFR-style synthetic graphs where the original FPA's
// k-core seed misses the true densest community. Runs baseline FPA,
// then expands 2-hop balls around top-degree vertices and runs MGPP
// locally, keeping the best.
AlgoResult runFPA3(const Graph& G, double tau);
