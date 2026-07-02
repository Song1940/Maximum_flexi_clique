# Maximum_Flexi_Clique

This is the implementation of the Maximum Flexi-Clique algorithms, which is described in the following paper submitted to SIGMOD 2027 (3rd round):
- Exact Maximum Flexi-Clique Computation under Non-Hereditary Feasibility and Connectivity Constraints

C++ implementation of exact algorithm (EBA) and heuristic (FPA) for finding the **maximum τ-flexi-clique** in an undirected graph.

A **τ-flexi-clique** of size *k* is a connected subgraph in which every node has degree at least ⌊k^τ⌋ within the subgraph (τ ∈ (0, 1)).

## Algorithms

| Binary | `--algo` | Description |
|--------|----------|-------------|
| `flexi_cum` | `ebac` | **EBA** — Efficient Branch-and-Bound (exact) |
| `flexi_cum` | `fpa` | **FPA** — Flexi-Prune Algorithm (heuristic) |
| `flexi_cum` | `npa` | **NPA** — Naive Peeling Algorithm (heuristic) Baseline algorithm |

The three algorithms share the `flexi_cum` binary and are selected via `--algo`. The build also produces ablation and baseline binaries (rule and branching-order variants, the SE/Naive baselines, and the query solver `qcs`).

## Requirements

- C++17
- CMake ≥ 3.14

## Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Usage

```
./build/flexi_cum --algo <npa|fpa|ebac> --file <path/to/network.dat> --tau <0.0–1.0> [--out <result.txt>]
```

**Example:**
```bash
./build/flexi_cum --algo ebac --file dataset/dolphin/network.dat --tau 0.9
```

**Output:**
```
=== Result ===
Algorithm : ebac
File      : dataset/dolphin/network.dat
Tau       : 0.9
Size      : 5
Time(s)   : 0.0004
Valid     : yes
Branches  : 0
```

## Input Format

Plain text edge list (whitespace-separated). Comment lines starting with `#` are ignored. Node IDs can be any non-negative integers; they are re-indexed to `[0, n)` internally.

```
# optional comment
0 1
0 2
1 3
...
```

## Datasets

Small graphs are included under `dataset/`, each as a `network.dat` edge list (with an optional `community.txt` for the recovery study). The five large graphs are not stored in the repository; download them from the sources below and convert them to the `u v` edge-list format:

| Dataset | Source |
|---------|--------|
| youtube | SNAP com-Youtube (https://snap.stanford.edu/data/com-Youtube.html) |
| livejournal | SNAP com-LiveJournal (https://snap.stanford.edu/data/com-LiveJournal.html) |
| pokec | SNAP soc-Pokec (https://snap.stanford.edu/data/soc-Pokec.html) |
| wiki | SNAP wiki-topcats (https://snap.stanford.edu/data/wiki-topcats.html) |
| florida | DIMACS 9th Challenge, Florida road network (http://www.diag.uniroma1.it/challenge9/download.shtml) |
