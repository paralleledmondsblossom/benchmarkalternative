#include <iostream>
#include <string>

enum class MaximumMatchingType {
    SEQUENTIAL_DFS,
    SEQUENTIAL_BFS,
    SEQUENTIAL_PF,
    SEQUENTIAL_PFP,
    SEQUENTIAL_HK,
    SEQUENTIAL_HK_DW,
    SEQUENTIAL_ABMP,
    SEQUENTIAL_ABMP_BFS,
    GPU_APFB_GPUBFS,
    GPU_APFB_GPUBFS_WR,
    GPU_APsB_GPUBFS,
    GPU_APsB_GPUBFS_WR
};

enum class InitialMatchingType {
    CHEAP_MATCHING,
    SK,
    SK_RAND,
    MIND_CHEAP,
    NO_INITIAL_MATCHING
};

const std::string MaximumMatchingTypeStrings[] = {
    "SequentialDFS",
    "SequentialBFS",
    "SequentialPF",
    "SequentialPFP",
    "SequentialHK",
    "SequentialHK_DW",
    "SequentialABMP",
    "SequentialABMP_BFS",
    "GPU-APFB_GPUBFS",
    "GPU-APFB_GPUBFS_WR",
    "GPU-APsB_GPUBFS",
    "GPU-APsB_GPUBFS_WR"
};

const std::string InitialMatchingTypeStrings[] = {
    "CheapMatching",
    "SK",
    "SK_rand",
    "mind_cheap",
    "no_initial_matching"
};
