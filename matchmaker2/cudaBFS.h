/*
---------------------------------------------------------------------
 This file is part of matchmaker framework
 Copyright (c) 2013,
 By:    Mehmet Deveci,
        Kamer Kaya,
        Bora Ucar,
        Umit V. Catalyurek
---------------------------------------------------------------------
 For license info, please see the LICENSE.txt file in
 the main directory.
---------------------------------------------------------------------
*/
//STARTLEVEL is set to -1

void swap_edges_GPUBFS_WR(
                int *_cmatch,
                int *_rmatch,
                int nc, 
                int nr, 
                int *_preced, 
                //int total_thread_num, 
                int blockD, 
                int threadDim, 
                int *bfs
                );

void GPUBFS_WR(
                  int level, 
                  int *_cxadj,
                  int *_cadj,
                  //int *cmatch,
                  int *_rmatch,
                  int nc, 
                  //int nr, 
                  int *_bfs, 
                  int *_preced, 
                  bool *_is_inserted,
                  bool *_non_matched_found, 
                  bool *_is_inserted2, 
                  //int total_thread_num, 
                  int *_root, 
                  int blockD,
                  int threadDim
                  );

void fixMatching(
                   int *_cmatch,
                   int *_rmatch,
                   int nr, int nc,
                   int *_bfs,
                   bool *_non_matched_found, 
                   //int total_thread_num,  
                   int blockD,
                   int threadDim
                   );

void init_BFSArray(
                     int nc,
                     int *_cmatch,
                     int *_bfs,
                     //int total_thread_num,
                     int blockD,
                     int threadDim 

                     );
                     
void initRoot_BFSArray(
                     int nc,
                     int *_cmatch,
                     int *_root, 
                     int *_bfs,
                     //int total_thread_num,  
                     int blockD,
                     int threadDim
                     );

void GPUBFS(int level, int *cxadj,
               int *cadj,
               int *rmatch,
               int nc,
               int *bfs,
               int *preced,
               bool *inserted,
               bool *non_matched_found, 
               bool *inserted2,
               int blockDim, 
               int threadDim
            );

void swap_edges_GPUBFS(
                int *cmatch,
                int *rmatch,
                int nc, 
                int nr, 
                int *preced,
                int blockDim, 
                int threadDim
                );

void fixMatching_initBFSArray(
                int *cmatch,
                int *rmatch,
                int nr, int nc,
                int *bfs, 
                bool *_non_matched_found,
                int blockDim, 
                int threadDim
                );
