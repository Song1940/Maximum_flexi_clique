#pragma once
#include "graph.h"
#include "npa.h"

// EBA: Efficient Branch-and-Bound Algorithm for Maximum Flexi-Clique.
//
// Pruning rules can be disabled at compile time via preprocessor defines:
//   -DNO_RULE1        disable Rule 1 (min adjusted degree check)
//   -DNO_RULE2        disable Rule 2 (upper-bound size check)
//   -DNO_RULE3        disable Rule 3 (diameter-based candidate pruning)
//   -DNO_RULE4        disable Rule 4 (degree-distance combined pruning)
//   -DNO_RULE5        disable Rule 5 (follower-based sibling pruning)
//   -DNO_RULE6        disable Rule 6 (global degree threshold pruning)
//   -DNO_DEGREE_ORDER disable degree-ascending branching order

AlgoResult runEBA(const Graph& G, double tau);

// EBA with cumulative CP-branching (exact maximum flexi-clique).
// Same maximum as runEBA / flexi_reference.py; branching differs (§ spec).
// Build with -DUSE_MARGIN_PIVOT to enable the margin-pivot ordering (default off).
AlgoResult runEBACumulative(const Graph& G, double tau);

// Seeded community search: maximum tau-flexi-clique containing `seed` within
// the given (local region) graph. Exact, anytime (best-so-far on timeout).
AlgoResult runEBASeeded(const Graph& region, int seed, double tau,
                        double timeout_sec);

// Exact branching baselines (R6-D3) — isolate the value of CP-branching.
// Same exact maximum, but Naive (include/exclude) / SE (set-enumeration) branching
// WITHOUT connectivity-by-construction, so connectivity is checked explicitly.
// AlgoResult.depth_prunes carries the number of explicit connectivity checks.
AlgoResult runNaiveBranching(const Graph& G, double tau);
AlgoResult runSEBranching(const Graph& G, double tau);
