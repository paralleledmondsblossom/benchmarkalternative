# alternativebenchmarks
This is a wrapper to call three matching solvers from one unified code.

Two initialization options are provided (ms-bfs-graft; matchmaker2) which perform bipartite matching with multiple threads or the GPU, respectively.  The first argument after the binary chooses one of these two.  The second argument determines whether to perform only initialization, only matching using third solver (Kececioglu) which performs matching using the Edmonds Blossom algorithm, or to use the bipartite matching as initialization to be completed by the Edmonds Blossom algorithm.  The third argument is the filename, and the next two arguments are initialization solver specific. 

Importantly, the bipartite codes don't produce simple matchings.  In the process of converting the graph to a bipartite graph, the same vertex appears in both sets.  This leads to the same vertex being matched more than one, which is invalid when the two sets are combined.  The initial matching size reported in the output file is a greedy packing of the matching created by the solver codes onto the original graph.

1. Buluc, Aydin, and Md Ariful Azad. Parallel Maximum Cardinality Matchings via Tree-Grafting. Lawrence Berkeley National Lab.(LBNL), Berkeley, CA (United States), 2016.
2. Deveci, Mehmet, et al. "GPU accelerated maximum cardinality matching algorithms for bipartite graphs." Euro-Par 2013 Parallel Processing: 19th International Conference, Aachen, Germany, August 26-30, 2013. Proceedings 19. Springer Berlin Heidelberg, 2013.
3. Kececioglu, John D., and A. Justin Pecqueur. "Computing Maximum-Cardinality Matchings in Sparse General Graphs." WAE. 1998.

## Installation

To get started with Alternative benchmarks for Maximum Cardinality Matching, follow these steps:

```bash
git clone https://github.com/paralleledmondsblossom/alternativebenchmarks.git
cd alternativebenchmarks

## Local Compilation
cd src
make

## Run
Usage: 
Wrapper args: a,b
a: ALGO [0: MS-BFS-GRAFT; 1: matchmaker2]
b: EXEC [0: JUST_DFS, 1: JUST_INIT, 2: FULL]
c: Filename (NOTE, this cannot be too long or buffer overflows)
MS-BFS args: d,e
d: Num Threads
e: KarpSipser Initialization [0: Serial, 1: Parallel]
matchmaker2 args:
d: Match type (variety of options)
e: Initial Match Type

MT: Maximum Matching Type
        0: Sequential DFS
        1: Sequential BFS
        2: Sequential PF
        3: Sequential PFP
        4: Sequential HK
        5: Sequential HK_DW
        6: Sequential ABMP
        7: Sequential ABMP_BFS
        8: GPU - APFB_GPUBFS
        9: GPU - APFB_GPUBFS_WR
        10: GPU - APsB_GPUBFS
        11: GPU - APsB_GPUBFS_WR
IMT: Initial Matching Type
        0: 
        0: Cheap Matching
        1: SK
        2: SK_rand
        3: mind_cheap
        >=4: no initial matching

binary         a b c                                       d  e
./src/matching 0 0 ../graphs/test_cases/luxembourg_osm.mtx 16 1
./src/matching 0 1 ../graphs/test_cases/luxembourg_osm.mtx 16 1
binary         a b c                                           d      e
./src/matching 1 0 IF=../graphs/test_cases/luxembourg_osm.mtx  MT=9   IMT=0
./src/matching 1 1 IF=../graphs/test_cases/luxembourg_osm.mtx  MT=9   IMT=0  

## Benchmarking
bash ./run_benchmarks.sh

## Input File

Matrix Market files are a standard format for representing sparse matrices. [Learn more](https://networkrepository.com/mtx-matrix-market-format.html).

## Output File

A CSV file is generated with the following columns:
- Solver used (arg1)
- Execution mode (arg2)
- Input File
- Number of Vertices
- Number of Edges
- Matching Size
- Parse Matrix Market File Time (milliseconds)
- Create Compressed Sparse Row Time (milliseconds)
- Initialization Time (seconds)
- Initialization Size
- CSR to Graph Time (seconds)
- SS DFS Time (Kececioglu) (seconds)
- Total Wall Clock Time (seconds)

The file is created with headers in the current working directory the executable is called from only once, and then new entries are appended to the existing file.

