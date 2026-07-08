#pragma once
#include "graph.h"
#include "npa.h"

// EBA: exact branch-and-bound algorithm for maximum flexi-clique.
//
// Pruning can be disabled at compile time for the ablation study:
//   -DNO_RULE1     disable Rule 1 (adjusted-degree pruning)
//   -DNO_DIAMETER  disable Rule 2 (diameter/distance pruning)
//   -DNO_RULE2     disable the size-bound termination check
//   -DNO_PIVOT     disable the width-based pivot


// EBA with cumulative CP-branching (exact maximum flexi-clique).
// Exact maximum flexi-clique (cumulative CP-branching).
// Build with -DUSE_MARGIN_PIVOT to enable the margin-pivot ordering (default off).
AlgoResult runEBACumulative(const Graph& G, double tau);

// Seeded community search: maximum tau-flexi-clique containing `seed` within
// the given (local region) graph. Exact, anytime (best-so-far on timeout).
AlgoResult runEBASeeded(const Graph& region, int seed, double tau,
                        double timeout_sec);

// Exact branching baselines — isolate the value of CP-branching.
// Same exact maximum, but Naive (include/exclude) / SE (set-enumeration) branching
// WITHOUT connectivity-by-construction, so connectivity is checked explicitly.
// AlgoResult.depth_prunes carries the number of explicit connectivity checks.
AlgoResult runNaiveBranching(const Graph& G, double tau);
AlgoResult runSEBranching(const Graph& G, double tau);
