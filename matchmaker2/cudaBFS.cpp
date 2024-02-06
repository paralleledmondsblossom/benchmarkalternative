#include "cudaBFS.h"
#include <iostream>
using namespace std;
__global__ void _GPUBFS(int level, 
                           int *cxadj,
                           int *cadj,
                           int *rmatch,
                           int nc, 
                           int *bfs, 
                           int *preced, 
                           bool *inserted,
                           bool *non_matched_found, 
                           bool *inserted2, 
                           int total_thread_num, 
                           int blockSize
                           ){
	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nc / total_thread_num;
  
  if (tx < nc % total_thread_num){
    process_cnt += 1;
  }
  bool inserted_any = false;
  bool nmf = false;
  
  int i = 0;
  for (i = 0; i < process_cnt; ++i){
    int myColumnVertex = i * total_thread_num + tx;
    
    if (bfs[myColumnVertex] == level){
      //      inserted_any = true;
      int j = cxadj[myColumnVertex];
      int end = cxadj[myColumnVertex + 1];
      for(; j < end; ++j){
        int neighborRow = cadj[j];
        int neighborColMatch = rmatch[neighborRow];
        
        if(neighborColMatch > -1){
          if(bfs[neighborColMatch] > -1){
            inserted_any = true;
            bfs[neighborColMatch] = level - 1;
            preced[neighborRow] = myColumnVertex;
          }
        } else if(neighborColMatch == -1){
          rmatch[neighborRow] = -2;
          preced[neighborRow] = myColumnVertex; 
          nmf = true;
        }
      }
    }
  }
  if (inserted_any){
    inserted[0] = true;
  }
  if(nmf){
    non_matched_found[0] = true;
  }
  if(tx == 0) {
    inserted2[0] = false;
  }
}


void GPUBFS(int level, int *cxadj,
               int *cadj,
               //int *cmatch,
               int *rmatch,
               int nc, 
               //int nr, 
               int *bfs, 
               int *preced, 
               bool *inserted,
               bool *non_matched_found, 
               bool *inserted2,
               int blockD, 
               int threadDim){

  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);

  _GPUBFS <<<dimGrid,dimBlock>>>
    (level, cxadj,
    cadj,
    rmatch,
    nc, 
    bfs, 
    preced, 
    inserted,
    non_matched_found, 
    inserted2, 
    blockD * threadDim,
    threadDim
   );
}




__global__ void _swap_edges_GPUBFS(
                           int *cmatch,
                           int *rmatch,
                           int nc, 
                           int nr, 
                           int *preced, 
                           int total_thread_num, 
                           int blockSize
                           ){
	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nr / total_thread_num;
  
  if (tx < nr % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  for (i = 0; i < process_cnt; ++i){
    int myRowVertex = i * total_thread_num + tx;
    if(rmatch[myRowVertex] == -2){
      int rowInd = myRowVertex;
      do{
        int matchedColumn = preced[rowInd];
        int mRow = cmatch[matchedColumn];
        if(preced[mRow] == matchedColumn){

          break;
        }
        
        cmatch[matchedColumn] = rowInd;
        rmatch[rowInd] = matchedColumn;
        rowInd = mRow;
      } while (rowInd != -1);
    }
  }
}

void swap_edges_GPUBFS(
                            int *cmatch,
                            int *rmatch,
                            int nc, 
                            int nr, 
                            int *preced,
                            int blockD, 
                            int threadDim){
  
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);

  _swap_edges_GPUBFS <<<dimGrid,dimBlock>>>(cmatch,
                                    rmatch,
                                    nc, 
                                    nr, 
                                    preced, blockD * threadDim,threadDim);
  
}



__global__ void _fixMatching_initBFSArray(
                            int *cmatch,
                            int *rmatch,
                            int nr, int nc,
                            int *bfs,
                            bool *_non_matched_found, 
                            int total_thread_num, 
                            int blockSize
                            ){

	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nr / total_thread_num;
  
  if (tx < nr % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  for (; i < process_cnt; ++i){
    int myRowVertex = i * total_thread_num + tx;
    int c = rmatch[myRowVertex];
    int r = cmatch[c];
    if (r != myRowVertex) 
      rmatch[myRowVertex] = -1;
  }
  
  process_cnt = nc / total_thread_num;
  
  if (tx < nc % total_thread_num){
    process_cnt += 1;
  }
  
  for (i = 0; i < process_cnt; ++i){
    int myColVertex = i * total_thread_num + tx;
    bfs[myColVertex] = cmatch[myColVertex];
  }
  if (tx == 0){
    _non_matched_found[0] = false;
  }
}

void fixMatching_initBFSArray(
                            int *cmatch,
                            int *rmatch,
                            int nr, int nc,
                            int *bfs,
                            bool *_non_matched_found,
                            int blockDim, 
                            int threadDim){
  
  
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockDim, 1,1);

  _fixMatching_initBFSArray <<<dimGrid,dimBlock>>>(cmatch,
                                     rmatch,
                                     nr, nc,
                                     bfs,
                                     _non_matched_found, 
                                     blockDim * threadDim, 
                                     threadDim);
  
}


__global__ void _initRoot_BFSArray(
                                 int nc,
                            int *cmatch,
                            int *root, int *bfs,
                            int total_thread_num, int blockSize
                            ){

	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nc / total_thread_num;
  
  if (tx < nc % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  for (; i < process_cnt; ++i){
    int myColVertex = i * total_thread_num + tx;
    if (cmatch[myColVertex] == -1){ 
      root[myColVertex] = myColVertex;
      bfs[myColVertex] = -1;
    }else{
      root[myColVertex] = -1;
      bfs[myColVertex] = 0;
    }
  }
}


void initRoot_BFSArray(
                      int nc,
                      int *_cmatch,
                      int *_root, int *_bfs,
                     //int total_thread_num,
                     int blockD,
                     int threadDim
                     ){
                     
  int total_thread_num = blockD * threadDim;
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  _initRoot_BFSArray<<<dimGrid,dimBlock>>>(
                                         nc,
                                         _cmatch,
                                         _root, _bfs,
                                         total_thread_num, threadDim);
}


__global__ void _init_BFSArray(
                            int nc,
                            int *cmatch,
                            int *bfs,
                            int total_thread_num, 
                            int blockSize
                            ){

	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nc / total_thread_num;
  
  if (tx < nc % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  for (i = 0; i < process_cnt; ++i){
    int myColVertex = i * total_thread_num + tx;
    bfs[myColVertex] = cmatch[myColVertex];
  }
}


void init_BFSArray(
                     int nc,
                     int *_cmatch,
                     int *_bfs,
                     //int total_thread_num,
                     int blockD,
                     int threadDim

                     ){
                     
  int total_thread_num = threadDim * blockD;
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  _init_BFSArray<<<dimGrid,dimBlock>>>(  nc,
                                         _cmatch,
                                         _bfs,
                                         total_thread_num, 
                                         threadDim);
}





/*
__global__ void _swap_edges_v2(
                            int *cmatch,
                            int *rmatch,
                            int nc, 
                            int nr, 
                            int *preced, int total_thread_num, int blockSize
                            ){
  //int total_thread_num = 1;
	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nr / total_thread_num;
  
  if (tx < nr % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  //for(; i < nr; ++i){
  for (i = 0; i < process_cnt; ++i){
    int myRowVertex = i * total_thread_num + tx;
    if(rmatch[myRowVertex] == -2){
      int rowInd = myRowVertex;
      do{
        int matchedColumn = preced[rowInd];
        int mRow = cmatch[matchedColumn];
        if(preced[mRow] == matchedColumn){
          break;
        }
        
        cmatch[matchedColumn] = rowInd;
        rmatch[rowInd] = matchedColumn;
        rowInd = mRow;
      } while (rowInd != -1);
    }
  }
}
*/
/*
void swap_edges_v2(
                    int *_cmatch,
                    int *_rmatch,
                    int nc, 
                    int nr, 
                   int *_preced, 
                   int threadDim, int blockD
                   ){
  
  int total_thread_num = blockD * threadDim;
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  
  _swap_edges_v2 <<<dimGrid,dimBlock>>>(_cmatch,
                                        _rmatch,
                                        nc, 
                                        nr, 
                                        _preced, total_thread_num,threadDim);
  
}
*/

__global__ void _fixMatching(
                            int *cmatch,
                            int *rmatch,
                            int nr, int nc,
                            int *bfs,
                            bool *_non_matched_found, 
                            int total_thread_num, 
                            int blockSize
                            ){

	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nr / total_thread_num;
  
  if (tx < nr % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  for (; i < process_cnt; ++i){
    int myRowVertex = i * total_thread_num + tx;
    int c = rmatch[myRowVertex];
    int r = cmatch[c];
    if (r != myRowVertex) 
      rmatch[myRowVertex] = -1;
  }
  
  if (tx == 0){
    _non_matched_found[0] = false;
  }
}

void fixMatching(
                    int *_cmatch,
                    int *_rmatch,
                    int nr, int nc,
                    int *_bfs,
                    bool *_non_matched_found, 
                    int blockD,
                    int threadDim
                    ){
  int total_thread_num = threadDim * blockD;
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  
  
  _fixMatching <<<dimGrid,dimBlock>>>(_cmatch,
                                        _rmatch,
                                        nr, nc,
                                        _bfs,_non_matched_found, 
                                        total_thread_num, 
                                        threadDim);
  
}

__global__ void _GPUBFS_WR(
                              int level, int *cxadj,
                              int *cadj,
                              //int *cmatch,
                              int *rmatch,
                              int nc, 
                              //int nr, 
                              int *bfs, 
                              int *preced, 
                              bool *inserted,
                              bool *non_matched_found, 
                              bool *inserted2, 
                              int total_thread_num, 
                              int blockSize,
                              int *root
                              ){
  
  
  
	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nc / total_thread_num;
  
  if (tx < nc % total_thread_num){
    process_cnt += 1;
  }


  bool inserted_any = false;
  bool nmf = false;
  int i = 0;
  for (i = 0; i < process_cnt; ++i){
    int myColumnVertex = i * total_thread_num + tx;
    //preced[0] = myColumnVertex - level;
    if (bfs[myColumnVertex] == level){
      
      
      int myRoot = root[myColumnVertex];
      if(bfs[myRoot ]> 0) continue;
      
      int j = cxadj[myColumnVertex];
      int end = cxadj[myColumnVertex + 1];
      for(; j < end; ++j){
        int neighborRow = cadj[j];
        int neighborColMatch = rmatch[neighborRow];
        
        if(neighborColMatch > -1){
          if(bfs[neighborColMatch] == 0){
            bfs[neighborColMatch] = level - 1;
            inserted_any = true;
            root[neighborColMatch] = myRoot;
            //todo may wanna check here.
            if (root[neighborColMatch] == myRoot)
              preced[neighborRow] = myColumnVertex;
          }
        } else if(neighborColMatch == -1){
          bfs[myRoot] = 1 +neighborRow;
          if(bfs[myRoot] == 1+neighborRow){
            
            //may want to check here.
            preced[neighborRow] = myColumnVertex; 
            rmatch[neighborRow] = -2;
            nmf = true;
          }
        }
      }
    }
  }
  
  //inserted[0] = false;
  if (inserted_any){
    inserted[0] = true;
  }
  if(nmf){
    non_matched_found[0] = true;
  }
  if(tx == 0) {
    inserted2[0] = false;
  }
   
}


void GPUBFS_WR(
                  int level, int *_cxadj,
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
                  ){
                  
  int total_thread_num = blockD * threadDim;
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  
  _GPUBFS_WR <<<dimGrid,dimBlock>>>
    (level, 
     _cxadj,
     _cadj,
     _rmatch,
     nc, 
     _bfs, 
     _preced, 
     _is_inserted,
     _non_matched_found, 
     _is_inserted2, 
     total_thread_num,
     threadDim, 
     _root
   );
}



__global__ void _swap_edges_GPUBFS_WR(
                               int *cmatch,
                               int *rmatch,
                               int nc, 
                               int nr, 
                               int *preced, 
                               int total_thread_num, 
                               int blockSize, 
                               int *bfs
                               ){
  //int total_thread_num = 1;
	int tx = blockIdx.x*blockSize + threadIdx.x;
  
  int process_cnt = nc / total_thread_num;
  
  if (tx < nc % total_thread_num){
    process_cnt += 1;
  }
  int i = 0;
  //for(; i < nr; ++i){
  for (i = 0; i < process_cnt; ++i){
    int myColVertex = i * total_thread_num + tx;
    //int myRowVertex = i * total_thread_num + tx;
    //int myRowVertex = i * total_thread_num + tx;
    //int myRowVertex = i;
    int rowInd = bfs[myColVertex] - 1;
    if(rowInd >= 0){
      do{
        int matchedColumn = preced[rowInd];
        int mRow = cmatch[matchedColumn];
        if(preced[mRow] == matchedColumn){
          
          break;
        }
        
        cmatch[matchedColumn] = rowInd;
        rmatch[rowInd] = matchedColumn;
        rowInd = mRow;
      } while (rowInd != -1);
    }
  }
}

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
                   ){
  int total_thread_num = blockD * threadDim;
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  
  _swap_edges_GPUBFS_WR <<<dimGrid,dimBlock>>>(_cmatch,
                                        _rmatch,
                                        nc, 
                                        nr, 
                                        _preced, 
                                        total_thread_num,
                                        threadDim,
                                        bfs
                                        );
  
}
