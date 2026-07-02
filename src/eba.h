#pragma once
#include "graph.h"
#include "npa.h"

// EBA: Efficient Branch-and-Bound Algorithm for Maximum Flexi-Clique.


AlgoResult runEBA(const Graph& G, double tau);

// EBA with cumulative CP-branching (exact maximum flexi-clique).
// Build with -DUSE_MARGIN_PIVOT to enable the margin-pivot ordering (default off).
AlgoResult runEBACumulative(const Graph& G, double tau);

// Seeded community search: maximum tau-flexi-clique containing `seed` within
// the given (local region) graph. Exact, anytime (best-so-far on timeout).
AlgoResult runEBASeeded(const Graph& region, int seed, double tau,
                        double timeout_sec);

// Same exact maximum, but Naive (include/exclude) / SE (set-enumeration) branching
// WITHOUT connectivity-by-construction, so connectivity is checked explicitly.
// AlgoResult.depth_prunes carries the number of explicit connectivity checks.
AlgoResult runNaiveBranching(const Graph& G, double tau);
AlgoResult runSEBranching(const Graph& G, double tau);
