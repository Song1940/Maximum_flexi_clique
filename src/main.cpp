#include "graph.h"
#include "npa.h"
#include "fpa.h"
#include "fpa3.h"
#include "eba.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <filesystem>
#include <unordered_set>
#include <algorithm>
#include <cstring>

namespace fs = std::filesystem;

static void printUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --algo <npa|fpa|eba>"
              << " --file <path/to/network.dat>"
              << " --tau <0.0-1.0>"
              << " [--out <output_file>]\n";
}

int main(int argc, char** argv) {
    std::string algo, file_path, out_file;
    double tau = 0.9;

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--algo") == 0 && i+1 < argc)
            algo = argv[++i];
        else if (std::strcmp(argv[i], "--file") == 0 && i+1 < argc)
            file_path = argv[++i];
        else if (std::strcmp(argv[i], "--tau") == 0 && i+1 < argc)
            tau = std::stod(argv[++i]);
        else if (std::strcmp(argv[i], "--out") == 0 && i+1 < argc)
            out_file = argv[++i];
    }

    if (algo.empty() || file_path.empty()) {
        printUsage(argv[0]);
        return 1;
    }

    // Load graph
    Graph G;
    try {
        G = Graph::loadFromFile(file_path);
    } catch (const std::exception& e) {
        std::cerr << "Error loading graph: " << e.what() << "\n";
        return 1;
    }

    // Run algorithm
    AlgoResult result;
    if (algo == "npa") {
        result = runNPA(G, tau);
    } else if (algo == "fpa") {
        result = runFPA(G, tau);
    } else if (algo == "fpa3") {
        result = runFPA3(G, tau);
    } else if (algo == "eba") {
        result = runEBA(G, tau);
    } else if (algo == "ebac") {
#ifdef EBA_CUM_AVAILABLE
        result = runEBACumulative(G, tau);
#else
        std::cerr << "Algorithm 'ebac' is not linked into this binary "
                     "(build the flexi_cum target).\n";
        return 1;
#endif
    } else if (algo == "se") {
#ifdef EBA_BASELINES_AVAILABLE
        result = runSEBranching(G, tau);
#else
        std::cerr << "Algorithm 'se' is not linked (build the flexi_baselines target).\n";
        return 1;
#endif
    } else if (algo == "naive") {
#ifdef EBA_BASELINES_AVAILABLE
        result = runNaiveBranching(G, tau);
#else
        std::cerr << "Algorithm 'naive' is not linked (build the flexi_baselines target).\n";
        return 1;
#endif
    } else {
        std::cerr << "Unknown algorithm: " << algo << "\n";
        printUsage(argv[0]);
        return 1;
    }

    // Validate result
    bool valid = false;
    if (!result.nodes.empty()) {
        valid = isFlexiClique(G, result.nodes, tau);
    }

    // Print to stdout
    std::cout << "=== Result ===\n"
              << "Algorithm : " << algo << "\n"
              << "File      : " << file_path << "\n"
              << "Tau       : " << tau << "\n"
              << "Size      : " << result.nodes.size() << "\n"
              << "Time(s)   : " << result.time_sec << "\n"
              << "Valid     : " << (valid ? "yes" : "no") << "\n";
    if (algo == "eba" || algo == "ebac" || algo == "se" || algo == "naive") {
        std::cout << "Branches  : " << result.branches << "\n";
        if (result.depth_prunes > 0) {
            const char* lbl = (algo == "se" || algo == "naive")
                                ? "ConnChecks" : "DepthPrune";
            std::cout << lbl << ": " << result.depth_prunes << "\n";
        }
    }

    // Write to output file if specified
    if (!out_file.empty()) {
        // Auto-create parent directory
        fs::path p(out_file);
        if (p.has_parent_path())
            fs::create_directories(p.parent_path());

        std::ofstream fout(out_file);
        if (!fout.is_open()) {
            std::cerr << "Cannot write to: " << out_file << "\n";
        } else {
            fout << "algorithm=" << algo << "\n"
                 << "tau=" << tau << "\n"
                 << "size=" << result.nodes.size() << "\n"
                 << "time=" << result.time_sec << "\n"
                 << "valid=" << (valid ? "1" : "0") << "\n"
                 << "branches=" << result.branches << "\n";

            // Write node list
            fout << "nodes=";
            for (int i = 0; i < (int)result.nodes.size(); ++i) {
                fout << result.nodes[i];
                if (i + 1 < (int)result.nodes.size()) fout << ",";
            }
            fout << "\n";
            std::cout << "Output written to: " << out_file << "\n";
        }
    }

    return 0;
}
