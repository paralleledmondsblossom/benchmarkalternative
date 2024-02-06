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

#include <iostream>
#include "MMArguments.h"
#include <sys/time.h>
#include "matchmaker.h"
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cudaBFS.h"
#include "MatchingTypeEnums.h"
#include <unistd.h>
#include "GreedyMatcher.cuh"



//#define profiling

typedef int IT;
typedef int WT;

#define MAXLINE 128*(1024*1024)
#define step 2.0f
#define STARTLEVEL -1
typedef struct pair {
	long long f;
	long long s;
} Pair;

template <class IT>
struct uSortItem
{
IT id;
//unsigned int val;
IT val;
};// uSortItem;


template <class IT>
void uqsort(IT n, uSortItem<IT> *arr);
void errexit(const char *fmt,...);
int * imalloc(long size, const char * msg);
void uprintf(const char *f_str,...);

void * umalloc(long size, const char * msg);



double rtclock();

void sort( int **pcxadj, int **pcadj,
                       int nc, int nedge
                       ){
  
  int *cxadj = *pcxadj; 
  int *cadj = *pcadj;
  
  uSortItem<int> *adj = new uSortItem<int>[nc];
  for(int i = 0; i < nc; ++i){
    adj[i].id = i;
    adj[i].val = cxadj[i + 1] - cxadj[i];
  }
  uqsort(nc, adj);
  
  
  cout << "nc:" << nc << endl;
  cout << "nedge:" <<  nedge << endl;
  
  int *ncxadj = new int [nc + 1];
  int *ncadj = new int [nedge ];
  
  int pp = 0;
  for(int i = 0; i < nc; ++i){
    int next = adj[i].id;
    ncxadj[i] = pp;
    for (int k = cxadj[next]; k < cxadj[next + 1]; ++k){
      ncadj[pp++] = cadj[k];
    }
  }
  ncxadj[nc] = pp;
  
  //cout << "pp:" << pp << endl;
  
  
  for(int i = 0; i < nc; ++i){
    int next = adj[i].id;
    int npincount = cxadj[next + 1 ] - cxadj[next];
    int ipincount = ncxadj[i+1] - ncxadj[i];
    /*
    if (npincount != ipincount){
      cout << "for i:" << i << " next:" << next<<" i:" << ipincount << " n:" << npincount << endl;
    }
    */
  }
  delete [] cxadj;
  delete []cadj;
  
  
   *pcxadj = ncxadj; 
   *pcadj = ncadj;
  
  delete[] adj;
}

void APFB_GPUBFS_WR( int *cxadj,
                       int *cadj,
                       int *cmatch,
                       int *rmatch,
                       
                       int nc, 
                       int nr, 
                       int *bfs, 
                       int *preced, 
                       bool *inserted,
                       bool *non_matched_found, bool *inserted2
                       ,int blockD, int threadDim,
                       int *root);

void APsB_GPUBFS_WR(int *cxadj,
                         int *cadj,
                         int *cmatch,
                         int *rmatch,
                         
                         int nc, 
                         int nr, 
                         int *bfs, 
                         int *preced, 
                         bool *_is_inserted,
                         bool *_non_matched_found, bool *_is_inserted2
                         ,int blockD, int threadDim,
                         int *root
                         );

void APsB_GPUBFS(int *_cxadj, int *_cadj, 
              int * _rmatch, 
              int *_cmatch,
              int nr, int nc, 
              int *_bfs, 
              int *_preced, 
              bool *_is_inserted, bool *_non_matched_found, 
              bool *_is_inserted2, int bd,int td);

void APFB_GPUBFS(int *_cxadj, int *_cadj, 
              int * _rmatch, 
              int *_cmatch,
              int nr, int nc, 
              int *_bfs, 
              int *_preced, 
              bool *_is_inserted, bool *_non_matched_found, 
              bool *is_inserted2, int bd, int td);

int pcmp(const void *v1, const void *v2);

void ReadGraphFromFile(FILE *fpin, 
                       IT *nrow, 
                       IT **pxadj, 
                       IT **padjncy, 
                       WT **padjncyw,
                       WT **ppvw, IT *nn);

template <typename IT>
void ReadGraphFromMMFile(
                         string filename, 
                         IT *nrow, 
                         IT *ncol,
                         IT *nnzero,
                         IT **pxadj, 
                         IT **padjncy);



void readGraph(MMArguments *mma, 
               IT *nv,
               IT *nc,
               IT *nn,
               IT **pxadj, 
               IT **padj);

void create_col_ptr(IT nv,
                    IT nc,
                    IT nn,
                    IT *pxadj, 
                    IT *padj,
                    IT **pcxadj,
                    IT **pcadj);

void checkSymmetry(IT nv, IT *pxadj,
                   IT *padj);

void initial_matching(IT nr,
                      IT nc,
                      IT *pxadj, 
                      IT *padj,
                      IT *pcxadj,
                      IT *pcadj,
                      IT *cmatch,
                      IT *rmatch,
                      int initial_matching_type);

IT match_count(IT nr,
               IT nc,
               IT *cmatch, IT *rmatch);

void maximal_matching(IT nr,
                      IT nc,
                      IT *pxadj, 
                      IT *padj,
                      IT *pcxadj,
                      IT *pcadj,
                      IT *cmatch,
                      IT *rmatch,
                      int matching_type,
                      IT *_rxadj, 
                      IT *_radj,
                      IT *_cxadj,
                      IT *_cadj,
                      IT *_cmatch,
                      IT *_rmatch,
                      bool *_is_inserted,
                      int *_bfs, 
                      int *_preced,
                      bool *_is_inserted2,
                      bool *_non_matched_found,
                      int blockDim, 
                      int threadDim,IT *_root_array);
#include <sys/time.h>
double getTimeOfDay2() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

extern "C"
int main_lib(int argc, char *argv[], FILE * fp, IT **cxadj, IT **cadj, IT **matching, IT*nr_ptr, IT*nc_ptr, IT*nn_ptr, int just_read_file, double *parse_graph_time, double *create_csr_time, double *init_time, int *init_match_count){
  MMArguments mma(argc, argv);
  
  IT nr, nc, nn;
  IT *rxadj;
  IT *radj;
  
  //IT *cxadj, *cadj;
  
  IT *_rxadj;
  IT *_radj;
  IT *_cxadj;
  IT *_cadj;
  IT *_cmatch;
  IT *_rmatch;
  bool *_is_inserted;
  IT *_bfs, *_preced;
  
  bool *_non_matched_found;
  bool *_is_inserted2;
  
  IT *_root_array = NULL;
  
  int match_types[11];
  /*
  match_types[3] = 1;
  match_types[4] = 1;
  match_types[5] = 1;
  */
  double start_parse_graph = getTimeOfDay2();
  readGraph(
            &mma, 
            nr_ptr, nc_ptr, nn_ptr,
            &rxadj, 
            &radj);
  double end_parse_graph = getTimeOfDay2();
  *parse_graph_time = end_parse_graph-start_parse_graph;
  nr = *nr_ptr;
  nc = *nc_ptr;
  nn = *nn_ptr;
  
  
  //cout << "creating columns:" << endl;
  double cbeg = rtclock();
  create_col_ptr(nr, nc, nn, rxadj, radj, cxadj, cadj);
  double cend = rtclock();
  cout << "col creation:" << cend - cbeg << endl;
  *create_csr_time = cend - cbeg;
  if (just_read_file)
    return 0;
  /*
  if (mma.sort){
    double sbeg = rtclock();
    sort(&cxadj, &cadj, nc, nn);
    double send = rtclock();
    cout << "Sorting:" << send - sbeg << endl;  
    
  }
  */
  {
#ifdef symmetry_check
  
  checkSymmetry(nr, rxadj, radj);
  for(IT i = 0; i < nc; ++i){
    if (cxadj[i] != rxadj[i]){
      cout << "i:" << i << " cxadj[i]:" << cxadj[i] << " rxadj[i]:" << rxadj[i] << endl;
      break;
    }
  }
  checkSymmetry(nc, cxadj,
                cadj);
#endif
  }
  *matching = new IT [nr];

  IT *cmatch = new IT[nc];
  IT *rmatch = new IT [nr];
  
  IT *cmatch_c = new IT[nc];
  IT *rmatch_c = new IT [nr];
  
  
  cout  << "nr:" << nr << endl
        << "nc:" << nc << endl
        << "nn:" << nn <<endl
        << "avg_edge:" << nn / float(nr) << endl;
  
  for (IT i = 0; i < nc; ++i) cmatch[i] = -1;
  for (IT i = 0; i < nr; ++i) rmatch[i] = -1;

  //cuda allocs
  size_t nnz_size = sizeof(IT) * nn;
  size_t nr_size = sizeof(IT) * (nr + 1);
  size_t nc_size = sizeof(IT) * (nc + 1);
  
  if(mma.match_type == 9 || mma.match_type == 8 || mma.match_type == -1 || 
     mma.match_type == 10 || mma.match_type == 11){
    if(!(
       //cudaMalloc((void**)&_rxadj, nr_size ) == cudaSuccess && 
       //cudaMalloc((void**)&_radj, nnz_size ) == cudaSuccess && 
       cudaMalloc((void**)&_cxadj, nc_size ) == cudaSuccess && 
       cudaMalloc((void**)&_cadj, nnz_size ) == cudaSuccess && 
       cudaMalloc((void**)&_cmatch, nc_size ) == cudaSuccess && 
       cudaMalloc((void**)&_rmatch, nr_size ) == cudaSuccess && 
       cudaMalloc((void**)&_is_inserted, sizeof(bool) ) == cudaSuccess&& 
       cudaMalloc((void**)&_is_inserted2, sizeof(bool) ) == cudaSuccess&& 
       cudaMalloc((void**)&_non_matched_found, sizeof(bool) ) == cudaSuccess&& 
       cudaMalloc((void**)&_bfs, nc_size ) == cudaSuccess && 
       cudaMalloc((void**)&_preced, nr_size ) == cudaSuccess && 
       cudaMalloc((void**)&_root_array, nc_size ) == cudaSuccess)
         //TODO -- sizeof this.
       ){
      printf("CANNOT ALLOCATE GPU MEMORY\n");
      exit(1);
    }
  }
  double imbegin = rtclock();
  initial_matching(nr, nc, rxadj, radj, *cxadj, *cadj, cmatch, rmatch, mma.initial_matching_type);
  double imend = rtclock();
  IT mc = match_count(nr,
                      nc,
                      cmatch, rmatch);
  
  int mRuns = 1;
  int rep = 1;
  if(mma.match_type == -1){
    match_types[0] = 8;
    match_types[1] = 8;
    match_types[2] = 9;
    match_types[3] = 9;
    match_types[4] = 10;
    match_types[5] = 10;
    
    match_types[6] = 11;
    match_types[7] = 11;
    
    match_types[8] = 4;
    match_types[9] = 3;
    mRuns = 10;
    if(mma.rep == -1){
      rep = 5;
    }
    else{
      rep = mma.rep;
    }
  } else{
    match_types[0] = mma.match_type;
    mRuns = 1;
    if(mma.rep == -1){
      rep = 1;
    }
    else{
      rep = mma.rep;
    }
  }
  for (int ii = 0; ii < mRuns; ++ii) {
    int mType = match_types[ii];
    for (int jj = 0; jj < rep; ++jj) {
      if(mType == 8 || mType == 9 || mType == 10 || mType == 11) {
        double bfsbeg = rtclock();
        cudaMemcpy(_cmatch, cmatch, nc_size, cudaMemcpyHostToDevice);
        cudaMemcpy(_rmatch, rmatch, nr_size, cudaMemcpyHostToDevice);
        //cudaMemcpy(_rxadj, rxadj, nr_size, cudaMemcpyHostToDevice);
        //cudaMemcpy(_radj, radj, nnz_size, cudaMemcpyHostToDevice);
        
        cudaMemcpy(_cxadj, *cxadj, nc_size, cudaMemcpyHostToDevice);
        cudaMemcpy(_cadj, *cadj, nnz_size, cudaMemcpyHostToDevice);
        
        cudaMemcpy(_bfs, cmatch, nc_size, cudaMemcpyHostToDevice);
        cudaMemcpy(_preced, rmatch, nr_size, cudaMemcpyHostToDevice);
        bool a[1];
        a[0] = false;
        
        cudaMemcpy(_non_matched_found, a, sizeof(bool), cudaMemcpyHostToDevice);
        cudaMemcpy(_is_inserted, a, sizeof(bool), cudaMemcpyHostToDevice);
        cudaMemcpy(_is_inserted2, a, sizeof(bool), cudaMemcpyHostToDevice);
        
        //if(mType == 9){
        //  cudaMemcpy(_root_array, cmatch, nc_size, cudaMemcpyHostToDevice);
        //}
        double bfsend = rtclock();
        cout << "memcopyTime:" << bfsend - bfsbeg << endl;
      }
      else{
        
        memcpy ( cmatch_c, cmatch, sizeof(IT) *nc );
        memcpy ( rmatch_c, rmatch, sizeof(IT) *nr );
      }
      
      int bd = mma.blockDim;
      int td = mma.threadDim;
      if (mma.match_type == -1 && ii % 2 == 0){
        td = 256;
        bd = nc / td;
        if (nc % td){
          bd++;
        }
      }
      if (mma.match_type == -1 && ii % 2 == 1){
        td = 256;
        bd = 256;
      }
      
      if (mma.match_type == 11 || mma.match_type == 8 || mma.match_type == 9 || mma.match_type == 10){
        if (bd == -1 && td == -1){
          td = 256;
          bd = nc / td;
          if (nc % td){
            bd++;
          }
        }
        if (bd == -1 ){
          bd = 256;
        }
        if (td == -1){
          td = 256;
        }
      }
      if (bd > 65535){
        bd = 65535;
      }
      double mbegin = rtclock();
      maximal_matching(nr, nc, rxadj, radj, 
                       *cxadj, *cadj, 
                       cmatch_c, rmatch_c, 
                       mType, 
                       _rxadj, _radj, 
                       _cxadj, _cadj, 
                       _cmatch, _rmatch, 
                       _is_inserted, 
                       _bfs, _preced,
                       _is_inserted2, 
                       _non_matched_found, bd, td,
                       _root_array);
      double mmend = rtclock();
      
      GreedyMatcher gm(nr,_cmatch,_rmatch);
      double mbegin2 = rtclock();
      int numAugmented = gm.maxMatch(*matching);
      double mmend2 = rtclock();
      if(mType == 8 || mType == 9 || mType == 10 ||  mType == 11) {
        cudaMemcpy(cmatch_c, _cmatch, sizeof(int) * nc, cudaMemcpyDeviceToHost);
        cudaMemcpy(rmatch_c,_rmatch, sizeof(int) * nr, cudaMemcpyDeviceToHost); 
      }
      
      IT mmc = match_count(nr,
                           nc,
                           cmatch_c, rmatch_c);
      cout << "Initial Match Time:" << imend - imbegin << endl;
      cout << "Initial Match Count:" << mc << endl;
      cout << "Maximal Match Time:" << mmend - mbegin << endl;
      cout << "Maximal Match Count:" << mmc << endl;
      cout << "Fixed Match Count:" << numAugmented << endl;
      *init_match_count=numAugmented;
      cout << "Fixed Match Time:" << mmend2 - mbegin2 << endl;
      *init_time=mmend2 - mbegin2;
      char inputFilename[500];
      char outputFilename[500];
      strcpy(outputFilename, "Results.csv");
      FILE *output_file;
      if (access(outputFilename, F_OK) == 0)
      {
          // file exists
          output_file = fopen(outputFilename, "a");
      }
      else
      {
          // file doesn't exist
          output_file = fopen(outputFilename, "w");
          fprintf(output_file, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "Filename","V","E","InitialMatchType", "MaxMatchType", "InitialMatchSize", "MaxMatchSize", "FixedMatchSize",  "InitialMatchTime", "MaxMatchTime", "FixedMatchTime");
      }
      if (argc>1){
          fprintf(output_file, "%s,%d,%d,%s,%s,%d,%d,%d,%f,%f,%f\n", mma.input_file.c_str(),nr,nn,
                                                                  InitialMatchingTypeStrings[mma.initial_matching_type].c_str(),
                                                                  MaximumMatchingTypeStrings[mType].c_str(),
                                                                  mc,
                                                                  mmc,
                                                                  numAugmented,
                                                                  imend - imbegin,
                                                                  mmend - mbegin,
                                                                  mmend2 - mbegin2);
      }
      fclose(output_file);

    }
  }
  
  
  
  delete []rxadj;
  delete []radj;
  
  //delete []cxadj;
  //delete []cadj;
  
  delete []cmatch;
  delete []rmatch;
  
  
  //cudaFree(_rxadj );
  //cudaFree(_radj );
  cudaFree(_cxadj );
  cudaFree(_cadj );
  cudaFree(_cmatch);
  cudaFree(_rmatch) ;
  return mma.match_type;
}



double rtclock()
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) cout << "Error return from gettimeofday:" << stat << endl;
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
//template <typename IT, WT>
void readGraph(MMArguments *mma, 
               IT *nv,
               IT *nc,
               IT *nn,
               IT **pxadj, 
               IT **padj){
  
  WT *padjncyw;
  WT *ppvw;FILE *fpin;
  switch(mma->input_type){
    case 0:
      cout << "Reading Metis Format:" << mma->input_file << endl;
      fpin = fopen(mma->input_file.c_str(), "r");
      if (fpin == NULL) {
        cout << "cannot read file" << endl;
        exit(1);
      }
      ReadGraphFromFile(fpin, 
                             nv, 
                             pxadj, 
                             padj, 
                             &padjncyw,
                             &ppvw, nn);
      *nc = *nv;
      cout << "nc:" << *nc << " nn:" << *nn << endl;
      
      break;
    case 1:
      cout << "Reading Matrix Market format graph:" << mma->input_file << endl;
      ReadGraphFromMMFile<IT>(mma->input_file, 
                              nv, nc, nn,
                              pxadj, 
                              padj);
      break;
    default:
      cout <<"Input type is not defined." << endl;
      exit(1);
  }
  
}

void create_col_ptr(IT nv,
                    IT nc,
                    IT nn,
                    IT *pxadj, 
                    IT *padj,
                    IT **pcxadj,
                    IT **pcadj){
  
  
  IT *cxadj = new IT[nc + 1];
  *pcxadj = cxadj;
  
  IT *cadj = new IT[nn];
  *pcadj = cadj;
  
  IT *p_it_nc_work = new IT[nc];
  
  //cout << "1" << endl;
  //cout << "nc:" << nc << " nv:" << nv << endl;
  memset(cxadj, 0, sizeof(IT) * (nc+1));
    //cout << "2" << endl;
    
  int a = 0;
  
  int prevbeg = 0;
  int prevend = 0;
  for (IT i = 0; i < nv; ++i){
    IT adb = pxadj[i];
    IT adbend = pxadj[i+1];
    
    //if (adbend > nn){
    //  
    //}
    
    for(IT j = adb; j < adbend; ++j){
      IT ac = padj[j];
      ++cxadj[ac];
      a++;
    }
    
  }
  
  int i = nv - 1;
    IT adb = pxadj[i];
    IT adbend = pxadj[i+1];
  //cout << "i:" << i << " asdb:" << adb << " adbend:" << adbend << " nn:" << nn << endl;
  
  //cout << "a:" << a << " nn:" << nn <<endl;
  
    //cout << "3" << endl;
  IT ps = 0;
  for (IT i = 0; i < nc; ++i){
    IT pinc = cxadj[i];
    cxadj[i] = ps;
    
    //if (ps > nn)
    //cout << "\ti:" << i << " cx:" << ps << endl;
    p_it_nc_work[i] = ps;
    ps += pinc;
  }
  cxadj[nc] = ps;
    //cout << "4" << endl;
  for (IT i = 0; i < nv; ++i){
    IT adb = pxadj[i];
    IT adbend = pxadj[i+1];
    //cout << "\t"<< i << " adb:" << adb << " abdend:"<< adbend<<endl;
    for(IT j = adb; j < adbend; ++j){
      //cout << "\t\tj:" << j << endl;
      IT ac = padj[j];
      //cout << "\t\tacs:" << ac << endl;
      IT adBegin  = p_it_nc_work[ac]++;
      
      //cout << "\t\tadBegin:" << adBegin << " nca:" << nn << endl;
      cadj[adBegin] = i;
      //cout << "\t\tset adBegin:" << i << endl;
    }
    
    //cout << "\t ends\t"<< i << endl;
  }
    //cout << "5" << endl;
#ifdef debug
  for (IT i = 1; i < nc + 1; ++i){
    if(p_it_nc_work[i - 1] != cxadj[i]){
      cout << "i:" << i <<  " " << p_it_nc_work[i - 1] << " " << cxadj[i] << endl;
      break;
    }
  }
#endif
  delete [] p_it_nc_work;
}

void checkSymmetry(IT nv, IT *pxadj,
                   IT *padj){
  for(IT i = 0; i < nv; ++i){
    IT adbegin = pxadj[i];
    IT adEnd = pxadj[i + 1];
    //cout << "i:" << i << " b:" << adbegin << " e:" << adEnd << endl;
    for (IT j = adbegin; j < adEnd; ++j){
      IT neighbor = padj[j];
      
      IT adbegin1 = pxadj[neighbor];
      IT adEnd1 = pxadj[neighbor+1];
      bool found = false;
      for (IT k = adbegin1; k < adEnd1; ++k){
        IT n1 = padj[k];
        if (n1 == i){
          found = true;
          break;
        }
      }
      if (!found){
        cout << "Graph is not symetric" << endl;
        return;
      }
    }
  }
  cout << "Graph is symetric" << endl;
}

void initial_matching(IT nr,
                      IT nc,
                      IT *pxadj, 
                      IT *padj,
                      IT *pcxadj,
                      IT *pcadj,
                      IT *cmatch,
                      IT *rmatch,
                      int initial_matching_type){

  switch (initial_matching_type) {
    case 0:
      old_cheap(pxadj,padj , rmatch, cmatch, nr, nc);
      //old_cheap(pcxadj,pcadj , cmatch, rmatch, nc, nr);
      break;
    case 1:
      sk_cheap(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;
    case 2:
      sk_cheap_rand(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;
    case 3:
      mind_cheap(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;
    case 4:
      break;
    default:
      break;
  }
}

IT match_count(IT nr,
               IT nc,
               IT *cmatch, IT *rmatch){
  
  IT mc = 0;
  for (IT i = 0; i < nc; ++i){
    if (cmatch[i] >= 0){
      IT r = cmatch[i];
      IT c = rmatch[r];
      if (c != i){
        cout << "NC Error at matching ";
        cout << "c:" << i << " mc:" << r << " mr:" << c << endl; 
        if(i > 100){
          break;
          exit(1);
        }
        exit(1);
      }
      ++mc;
    }
  }
  
  IT mr = 0;
  for (IT i = 0; i < nr; ++i){
    if (rmatch[i] >= 0){
      IT c = rmatch[i];
      IT r = cmatch[c];
      if (r != i){
        cout << "NR Error at matching ";
        cout << "r:" << i << " mr:" << c << " mc:" << r << endl;
        if(i > 100){
        break;
        exit(1);
        }
        exit(1);
      }
      ++mr;
    }
  }
  if (mr != mc){
    cout << "mr:" << mr << " mc:" << mc << endl;
  }
  return mc;
}


void maximal_matching(IT nr,
                      IT nc,
                      IT *pxadj, 
                      IT *padj,
                      IT *pcxadj,
                      IT *pcadj,
                      IT *cmatch,
                      IT *rmatch,
                      int matching_type,
                      IT *_rxadj, 
                      IT *_radj,
                      IT *_cxadj,
                      IT *_cadj,
                      IT *_cmatch,
                      IT *_rmatch,
                      bool *_is_inserted,
                      int *_bfs, 
                      int *_preced,
                      bool *_is_inserted2,
                      bool *_non_matched_found,
                      int bd, 
                      int td,IT *_root_array){
  

  switch (matching_type) {
    case 0:
      cout  << endl << endl << "###############################" 
      << endl << "Running DFS" << endl;
      match_dfs(pcxadj,pcadj , cmatch, rmatch, nc, nr);
      break;
    case 1:
      cout  << endl << endl << "###############################"  
      << endl << "Running BFS" << endl;
      match_bfs(pcxadj,pcadj , cmatch, rmatch, nc, nr);
      break;
    case 2:
      cout  << endl << endl << "###############################"  
      << endl << "Running PF" << endl;
      match_pf(pcxadj,pcadj , cmatch, rmatch, nc, nr);
      break;
    case 3:
      cout  << endl << endl << "###############################" 
      << endl << "Running PFP" << endl;
      match_pf_fair(pcxadj,pcadj , cmatch, rmatch, nc, nr);
      
      break;
    case 4:
      
      cout  << endl << endl << "###############################" 
      << endl << "Running HK" << endl;
      match_hk(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;    
    case 5:
      cout  << endl << endl << "###############################" 
      << endl << "Running HK_DW" << endl;
      match_hk_dw(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;
    case 6:
      cout  << endl << endl << "###############################" 
      << endl << "Running ABMP" << endl;
      match_abmp(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;
    case 7:
      cout  << endl << endl << "###############################" 
      << endl << "Running ABMP_BFS" << endl;
      match_abmp_bfs(pcxadj,pcadj , pxadj, padj, cmatch, rmatch, nc, nr);
      break;
    case 8:
      cout  << endl << endl << "###############################" 
      << endl << "Running APFB_GPUBFS with " << bd <<" blocks and "<< td << " threads" << endl;      
      APFB_GPUBFS(_cxadj, _cadj, 
                    _rmatch, 
                    _cmatch,
                    nr, nc, 
                    _bfs, 
                    _preced, 
                    _is_inserted, _non_matched_found, 
                    _is_inserted2, bd, td);
      break;
    case 9:
      cout  << endl << endl << "###############################" 
      << endl << "Running APFB_GPUBFS_WR with " << bd <<" blocks and "<< td << " threads" << endl;
      APFB_GPUBFS_WR(  
                        _cxadj,
                        _cadj,
                        _cmatch,
                        _rmatch,
                        
                        nc, 
                         nr, 
                        _bfs, 
                        _preced, 
                        _is_inserted,
                        _non_matched_found, 
                        _is_inserted2
                        ,bd, td,
                        _root_array);
      break;
    case 10:
      cout  << endl << endl << "###############################" 
      << endl << "Running APsB_GPUBFS (early exit) with " << bd <<" blocks and "<< td << " threads" << endl;      
      APsB_GPUBFS(_cxadj, _cadj, 
               _rmatch, 
               _cmatch,
               nr, nc, 
               _bfs, 
               _preced, 
               _is_inserted, _non_matched_found, 
               _is_inserted2, bd, td);
      break;
    case 11:
      cout  << endl << endl << "###############################" 
      << endl << "Running APsB_GPUBFS_WR with " << bd <<" blocks and "<< td << " threads" << endl;
      APsB_GPUBFS_WR(  
                        _cxadj,
                        _cadj,
                        _cmatch,
                        _rmatch,
                        
                        nc, 
                        nr, 
                        _bfs, 
                        _preced, 
                        _is_inserted,
                        _non_matched_found, _is_inserted2
                        ,bd, td,
                        _root_array
                        );

      break;
      
    default:
      break;
  }
  
}

int pcmp(const void *v1, const void *v2){
	long long diff = (((Pair *)v1)->f - ((Pair *)v2)->f);
	if (diff != 0)
		return diff;
	else
		return (((Pair *)v1)->s - ((Pair *)v2)->s);
}


/* reads the Matrix Market format graph */

template <typename IT>
void ReadGraphFromMMFile(
                         string filename, 
                         IT *nrow, 
                         IT *ncol,
                         IT *nnzero,
                         IT **pxadj, 
                         IT **padjncy)
{
  
	
	string line;
  bool is_symmetric = false;
  ifstream filestr (filename.c_str(), fstream::in);
  if (!filestr.is_open()){
    cout << "File " << filename << " cannot be opened." << endl;
    exit(1);
  }
  do{
    getline(filestr, line);
    if (line[0] == '%' && line[1] == '%'){
      
      stringstream ss (stringstream::in | stringstream::out);
      ss << line;
      while(!ss.eof()){
        string el;
        ss >> el;
        if (el == "symmetric"){
          cout << "Symetric graph" << endl;
          is_symmetric = true;
          break;
        }
      }
    }
    if(filestr.eof()){
      cout << "unexpected EOF" << endl; exit(1);
    }
  }while(line[0] == '%');
  
  
  IT nr, nc, nnz; 
  stringstream ss (stringstream::in | stringstream::out);
  ss << line;
  ss >> nr >> nc >> nnz;
  
#ifdef debug
  cout << "nr:" << nr << " nc:" << nc << " nnz:" << nnz << endl;
#endif 
  
  int scale = 1;
  if (is_symmetric) scale = 2;
  
  nnz = nnz * scale;
  Pair *coords = (Pair*) malloc(sizeof(Pair) * nnz );
  
  IT i = 0;
  while(1){
    getline(filestr, line);
    if(filestr.eof()){
      break;
    }  
    IT v1, v2;
    stringstream ss (stringstream::in | stringstream::out);
    ss << line;
    ss >> v1 >> v2 ;
    
    
    coords[i].f = v1 - 1;
    coords[i].s = v2 - 1;
    ++i;
    if(is_symmetric){
      
      coords[i].f = v2 - 1;
      coords[i].s = v1 - 1;
      ++i;    
    }
  }
  filestr.close();
	qsort(coords, nnz, sizeof(Pair), pcmp);
  
  
	IT onnz = 1;
	for(IT i = 1; i < nnz; i++) {
		if(coords[i].f != coords[onnz-1].f || coords[i].s != coords[onnz-1].s) {
      //cout << "Eliminates" << endl;
			coords[onnz].f = coords[i].f;
			coords[onnz++].s = coords[i].s;
		}
	}
  
  nnz = onnz;
	*nrow = nr;
  *ncol = nc;
  *nnzero = nnz;
	
  IT *xadj = *pxadj = (int*) malloc((nr+1) * sizeof(IT));
	IT *adj = *padjncy = (int*) malloc(onnz * sizeof(IT));
  
  IT prevSource = 0;
  xadj[0] = 0;
  IT adjInd = 0;
  for (IT i = 0; i < nnz; ++i){
    IT target = coords[i].s;
    IT source = coords[i].f;
    
    while(prevSource < source){
      //cout << "xadj:" << prevSource << " " << adjInd << endl;
      xadj[++prevSource] = adjInd;
/*
      if (prevSource > nr - 5){
        cout << "xadj[prevSource]:" << xadj[prevSource] << endl;
        cout << "source:" << source << " prevsource:" << prevSource << endl;
    }
    */
    }
    adj[adjInd++] = target;
  }
  
  while(prevSource < nr){
    xadj[++prevSource] = adjInd;
  }
    
  xadj[nr] = adjInd;
  
  /*
  cout << "nr:" << nr << endl;
  cout << "xasdj:" << xadj[nr-1] << " xadj+1:" << xadj[nr] << endl;
  cout << "xasdj:" << xadj[nr-3] << " xadj-1:" << xadj[nr-2] << endl;
  */
  if (adjInd != nnz){
    cout << "error: adjInd:" << adjInd << " nnz:" << nnz << endl; 
    exit(1);
  }
	free(coords);
  
	return;
}



void ReadGraphFromFile(FILE *fpin, IT *numofvertex, IT **pxadj, IT **padjncy, WT **padjncyw,
                       WT **ppvw, IT *nn) {
  
	IT *xadj, *adjncy,  nvtxs, nedges, fmt, readew, readvw, edge, i, k, ncon;
	WT*adjncyw=NULL, *pvw=NULL;
	char *line;
  
	line = (char *)malloc(sizeof(char)*(MAXLINE+1));
  
	do {
		fgets(line, MAXLINE, fpin);
	} while (line[0] == '%' && !feof(fpin));
  
	if (feof(fpin))
		errexit("empty graph!!!");
  
	fmt = 0;
	{
		std::string s = line;
		std::stringstream ss (s);
		ss>>nvtxs>>nedges>>fmt>>ncon;
	}
	*numofvertex = nvtxs;
  
	readew = (fmt%10 > 0);
	readvw = ((fmt/10)%10 > 0);
	if (fmt >= 100)
		errexit("invalid format");
  
	nedges *=2;
  
  *nn = nedges;
	xadj = *pxadj = imalloc(nvtxs+2, "ReadGraph: xadj");
	adjncy = *padjncy = imalloc(nedges, "ReadGraph: adjncy");
	if (padjncyw)
		adjncyw = *padjncyw = imalloc(nedges, "ReadGraph: adjncyw");
	if (ppvw)
		pvw = *ppvw = imalloc(nvtxs+1, "ReadGraph: adjncyw");
  
	for (xadj[0]=0, k=0, i=0; i<nvtxs; i++) {
		char *oldstr=line, *newstr=NULL;
		int  ewgt=1, vw=1;
    
		do {
			fgets(line, MAXLINE, fpin);
		} while (line[0] == '%' && !feof(fpin));
    
		if (strlen(line) >= MAXLINE-5)
			errexit("\nBuffer for fgets not big enough!\n");
    
		if (readvw) {
			vw = (int)strtol(oldstr, &newstr, 10);
			oldstr = newstr;
		}
    
		if (ppvw)
			pvw[i] = vw;
    
		for (;;) {
			edge = (int)strtol(oldstr, &newstr, 10) -1;
			oldstr = newstr;
      
			if (readew) {
				ewgt = (int)strtol(oldstr, &newstr, 10);
				oldstr = newstr;
			}
      
			if (edge < 0) {
				break;
			}
      
			if (edge==i)
				errexit("Self loop in the graph for vertex %d\n", i);
      
      
			bool flag = false;
			for (int j = xadj[i]; j < k; j++) {
				if (adjncy[j] == edge) {
					flag = true;
					break;
				}
        
			}
      
			if (!flag) {
				adjncy[k] = edge;
        
				if (padjncyw)
					adjncyw[k] = ewgt;
				k++;
			}
		}
		xadj[i+1] = k;
	}
  
  //	if (k != nedges)
  //		errexit("k(%d)!=nedges(%d) and i:%d", k, nedges,i);
  
	free(line);
  
	return;
}


void APFB_GPUBFS(int *_cxadj, int *_cadj, 
              int * _rmatch, 
              int *_cmatch,
              int nr, int nc, 
              int *_bfs, 
              int *_preced, 
              bool *_is_inserted, bool *_non_matched_found, 
              bool *_is_inserted2, int bd,int td){
  bool is_inserted[1];
  *is_inserted = true;
  
  bool  non_matched_found[1];
  int level = STARTLEVEL;

  init_BFSArray(
                     nc,
                     _cmatch,
                     _bfs,
                     //total_thread_num,
                     bd,
                     td
                     );
#ifdef profiling
  double bfs_ = 0, swap = 0, init = 0;
  int totalbfscall = 0;
#endif
  while(1){
#ifdef profiling
    double bfsbeg = rtclock();
#endif
    while(1){
      
      GPUBFS (  level, 
                 _cxadj,
                 _cadj,
                 _rmatch,
                 nc, 
                 _bfs, 
                 _preced, 
                 _is_inserted,
                 _non_matched_found, 
                 _is_inserted2,
                 bd,
                 td
                 );
      
      
      cudaMemcpy(is_inserted,_is_inserted, sizeof(bool), cudaMemcpyDeviceToHost);
      if(!is_inserted[0]){
        cudaMemcpy(non_matched_found,_non_matched_found, sizeof(bool), cudaMemcpyDeviceToHost);
        break;
      }
      bool * tmp = _is_inserted;
      _is_inserted = _is_inserted2;
      _is_inserted2 = tmp;
      level--;
    }
#ifdef profiling
    double bfsend = rtclock();
    bfs_ += (bfsend - bfsbeg);
#endif 
    if(!non_matched_found[0]){
      break;
    }
    
    swap_edges_GPUBFS(
               _cmatch,
               _rmatch,
               nc, 
               nr, 
               _preced,bd,td
               );
    //double swapend = rtclock();              
    
    fixMatching_initBFSArray(
               _cmatch,
               _rmatch,
               nr,  nc,
               _bfs,
               _non_matched_found,bd,td
               );
    //double initend = rtclock();              
    
    //swap += (swapend - bfsend);
    //init += (initend - swapend);
#ifdef profiling
    cout << "\tbfs level:" << level << endl;
//    cout << "\tbfs time:" << (bfsend - bfsbeg) << endl;
    totalbfscall -= level;
#endif
    level = STARTLEVEL;
  }
  
#ifdef profiling
  cout << "bfs:" << bfs_ << endl;  
  cout << "bfscount:" << totalbfscall << endl;
#endif  

  
  
}



void APsB_GPUBFS(int *_cxadj, int *_cadj, 
              int * _rmatch, 
              int *_cmatch,
              int nr, int nc, 
              int *_bfs, 
              int *_preced, 
              bool *_is_inserted, bool *_non_matched_found, 
              bool *_is_inserted2, int bd,int td){
  bool is_inserted[1];
  *is_inserted = true;
  
  bool  non_matched_found[1];
  int level = STARTLEVEL;
#ifdef profiling
  double bfs_ = 0, swap = 0, init = 0;
  int totalbfscall = 0;
#endif

  bool done = false;
  while(1){
#ifdef profiling
    double bfsbeg = rtclock();
#endif


  init_BFSArray(
                     nc,
                     _cmatch,
                     _bfs,
                     //total_thread_num,
                     bd,
                     td 
                     );
    while(1){
      
      GPUBFS (level, 
                 _cxadj,
                 _cadj,
                 _rmatch,
                 nc, 
                 _bfs, 
                 _preced, 
                 _is_inserted,_non_matched_found, 
                 _is_inserted2,bd,td);
      
      
      cudaMemcpy(is_inserted,_is_inserted, sizeof(bool), cudaMemcpyDeviceToHost);
      cudaMemcpy(non_matched_found,_non_matched_found, sizeof(bool), cudaMemcpyDeviceToHost);
      if(!is_inserted[0]  && !non_matched_found[0]){
        done = true;
        break;
      }
        
      if(!is_inserted[0]  || non_matched_found[0]){
        break;
      }
      
      bool * tmp = _is_inserted;
      _is_inserted = _is_inserted2;
      _is_inserted2 = tmp;
      level--;
    }
#ifdef profiling
    double bfsend = rtclock();
    bfs_ += (bfsend - bfsbeg);
#endif
    if(done){
      break;
    }
    
    swap_edges_GPUBFS(
               _cmatch,
               _rmatch,
               nc, 
               nr, 
               _preced,bd,td
               );
    
    
    fixMatching_initBFSArray(
               _cmatch,
               _rmatch,
               nr,  nc,
               _bfs,
               _non_matched_found,bd,td
               );
    
    
#ifdef profiling
    cout << "\tbfs level:" << level << endl;
//    cout << "\tbfs time:" << (bfsend - bfsbeg) << endl;
    totalbfscall -= level;
#endif
    level = STARTLEVEL;
  }
#ifdef profiling
  cout << "bfs:" << bfs_ << endl;  
  cout << "bfscount:" << totalbfscall << endl;
#endif  

  
}


void APFB_GPUBFS_WR(int *cxadj,
                       int *cadj,
                       int *cmatch,
                       int *rmatch,
                       
                       int nc, 
                       int nr, 
                       int *bfs, 
                       int *preced, 
                       bool *_is_inserted,
                       bool *_non_matched_found, bool *_is_inserted2
                       ,int blockD, int threadDim,
                       int *root){
  
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  
  bool is_inserted[1];
  bool  non_matched_found[1];
  int level = STARTLEVEL;
  
  int total_thread_num = threadDim * blockD;

#ifdef profiling
  double bfs_ = 0, swap = 0, init = 0;
  int totalbfscall = 0;
#endif

  while(1){
    initRoot_BFSArray(
                    nc,
                    cmatch,
                    root, bfs,
                    //total_thread_num, 
                    blockD,
                    threadDim
                    );
    
#ifdef profiling
    double bfsbeg = rtclock();
#endif
    while(1){
      
      GPUBFS_WR(
                  level, cxadj,
                  cadj,
                  rmatch,
                  nc, 
                  bfs, 
                  preced, 
                  _is_inserted,
                  _non_matched_found, 
                  _is_inserted2, 
                  //total_thread_num, 
                  root,
                  blockD,
                  threadDim
                );
      
      
      cudaMemcpy(is_inserted,_is_inserted, sizeof(bool), cudaMemcpyDeviceToHost);
      if(!is_inserted[0]){
        cudaMemcpy(non_matched_found,_non_matched_found, sizeof(bool), cudaMemcpyDeviceToHost);
        break;
      }
      
      bool * tmp = _is_inserted;
      _is_inserted = _is_inserted2;
      _is_inserted2 = tmp;
      level--;
    }

#ifdef profiling
    double bfsend = rtclock();
    bfs_ += (bfsend - bfsbeg);
#endif
    
    if(!non_matched_found[0]){
      break;
    }
    
            
    
    swap_edges_GPUBFS (cmatch,
                   rmatch,
                   nc,
                   nr,
                   preced, 
                   blockD,
                   threadDim
                   );
    
    /*
    swap_edges_GPUBFS_WR(cmatch,
                    rmatch,
                    nc,
                    nr,
                    preced, 
                    //total_thread_num,
                    blockD, 
                    threadDim,
                    bfs);
    */
    fixMatching (cmatch,
                  rmatch,
                  nr, nc,
                  bfs,
                  _non_matched_found, 
                  //total_thread_num, 
                  blockD,
                  threadDim
                  );




#ifdef profiling
    cout << "\tbfs level:" << level << endl;
    totalbfscall -= level;
#endif
    level = STARTLEVEL;
  }
#ifdef profiling
  cout << "bfs:" << bfs_ << endl;  
  cout << "bfscount:" << totalbfscall << endl;
#endif  

  
}



void APsB_GPUBFS_WR(int *cxadj,
                       int *cadj,
                       int *cmatch,
                       int *rmatch,
                       int nc, 
                       int nr, 
                       int *bfs, 
                       int *preced, 
                       bool *_is_inserted,
                       bool *_non_matched_found, 
                       bool *_is_inserted2,
                       int blockD, 
                       int threadDim,
                       int *root
                       ){
  
  dim3 dimBlock(threadDim,1,1);
  dim3 dimGrid(blockD, 1,1);
  
  bool is_inserted[1];
  bool  non_matched_found[1];
  int level = STARTLEVEL;
  
  int total_thread_num = threadDim * blockD;
  bool done = false;
  
  
#ifdef profiling
  double bfs_ = 0, swap = 0, init = 0;
  int totalbfscall = 0;
#endif
  while(1){
    initRoot_BFSArray(
                    nc,
                    cmatch,
                    root, bfs,
                    //total_thread_num, 
                    blockD,
                    threadDim);
    
   
   
#ifdef profiling
    double bfsbeg = rtclock();
#endif
    while(1){
      
      GPUBFS_WR(
                   level, cxadj,
                   cadj,
                   rmatch,
                   nc, 
                   bfs, 
                   preced, 
                   _is_inserted
                   ,_non_matched_found, 
                   _is_inserted2, 
                   //total_thread_num, 
                   root,
                   blockD,
                   threadDim
                   );
      
     
      cudaMemcpy(is_inserted,_is_inserted, sizeof(bool), cudaMemcpyDeviceToHost);
      cudaMemcpy(non_matched_found,_non_matched_found, sizeof(bool), cudaMemcpyDeviceToHost);
      if(!is_inserted[0]  && !non_matched_found[0]){
        done = true;
        break;
      }
      if(!is_inserted[0]  || non_matched_found[0]){
        break;
      }
      
      bool * tmp = _is_inserted;
      _is_inserted = _is_inserted2;
      _is_inserted2 = tmp;
      level--;
  }
    
    
#ifdef profiling
    double bfsend = rtclock();
    bfs_ += (bfsend - bfsbeg);
#endif
    
    if(done){
     
      break;
    }
    
  
    
    /*
    swap_edges_GPUBFS (cmatch,
                   rmatch,
                   nc,
                   nr,
                   preced, 
                   blockD,
                   threadDim);
    */
    swap_edges_GPUBFS_WR(cmatch,
                    rmatch,
                    nc,
                    nr,
                    preced, 
                    //total_thread_num,
                    blockD, 
                    threadDim,
                    bfs
                    );
                    
                   
    
    fixMatching (cmatch,
                   rmatch,
                   nr, nc,
                   bfs,
                   _non_matched_found, 
                   //total_thread_num, 
                   blockD,
                   threadDim);


#ifdef profiling
    cout << "\tbfs level:" << level << endl;
    totalbfscall -= level;
#endif
    level = STARTLEVEL;
  }
#ifdef profiling
  cout << "bfs:" << bfs_ << endl;  
  cout << "bfscount:" << totalbfscall << endl;
#endif  

}



#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

#define U_MB (1024*1024)
#define M             7
#define NSTACK        50

template <class IT>
void uqsort(IT n, uSortItem<IT> * arr)
{
IT         i, ir=n, j, k, l=1;
IT         jstack=0, istack[NSTACK], aval;
uSortItem<IT>    a, temp;

--arr;
for (;;)
    {
    if (ir-l < M)
	{
	for (j=l+1;j<=ir;j++)
	    {
	    a=arr[j];
	    aval = a.val;
	    for (i=j-1;i>=1;i--)
		{
		if (arr[i].val <= aval)
		    break;
		arr[i+1] = arr[i];
		}
	    arr[i+1]=a;
	    }
	if (jstack == 0)
	    break;
	ir=istack[jstack--];
	l=istack[jstack--];
	}
    else
	{
	k=(l+ir) >> 1;
	SWAP(arr[k],arr[l+1])
	  if (arr[l+1].val > arr[ir].val)
	      {
	      SWAP(arr[l+1],arr[ir])
		}
	if (arr[l].val > arr[ir].val)
	    {
	    SWAP(arr[l],arr[ir])
	      }
	if (arr[l+1].val > arr[l].val)
	    {
	    SWAP(arr[l+1],arr[l])
	      }
	i=l+1;
	j=ir;
	a=arr[l];
	aval = a.val;
	for (;;)
	    {
	    do i++; while (arr[i].val < aval);
	    do j--; while (arr[j].val > aval);
	    if (j < i) break;
	    SWAP(arr[i],arr[j]);
	    }
	arr[l]=arr[j];
	arr[j]=a;
	jstack += 2;
	if (jstack > NSTACK){
		cout << "uqsort: NSTACK too small in sort." << endl;
		exit(1);
	}
	if (ir-i+1 >= j-l)
	    {
	    istack[jstack]=ir;
	    istack[jstack-1]=i;
	    ir=j-1;
	    }
	else
	    {
	    istack[jstack]=j-1;
	    istack[jstack-1]=l;
	    l=i;
	    }
	}
    }
}


void uprintf(const char *f_str,...)
{
va_list argp;
#ifdef _DEBUG
 char *name;
 static int cnt=0;
 FILE *ofp;
#endif

fflush(stdout);
fflush(stderr);
va_start(argp, f_str);
vfprintf(stdout, f_str, argp);
va_end(argp);
fflush(stdout);

/* ugly: repeating udbgprintf here because I don't know how to call that */
#ifdef _DEBUG
/* --------- print to debug file --------- */
 if (*_dbg_file_name)
     name = _dbg_file_name;
 else
     name = DEFAULT_DEBUG_FILE;
 
 if (cnt++==0)
    {
    time_t curtime;
    char    tst[26];

    ofp = ufopen(name, "w", "dbgprint");
    time(&curtime);
    strcpy(tst, ctime(&curtime));
    tst[24] = 0;
    fprintf(ofp, "\n\n\n----------------- %s -----------------\n\n", tst);
    }
 else
     ofp = ufopen(name, "a", "dbgprint");
va_start(argp, f_str);
vfprintf(ofp, f_str, argp);
va_end(argp);
ufclose(ofp);
#endif
}

void errexit(const char * f_str,...)
{
va_list argp;

fflush(stdout);
fflush(stderr);
fprintf(stderr, "\n****** Error:\n");
va_start(argp, f_str);
vfprintf(stderr, f_str, argp);
va_end(argp);

fprintf(stderr," ******\n");
fflush(stderr);

uprintf("Error in Execution\n");

exit(1);
}


void * umalloc(long size, const char * msg)
{
void * ptr;
static long __u_tot_size=0;

if (size == 0)
    return NULL;

ptr = (void *) malloc(size);
if (!ptr)
    {
    errexit("Memory allocation failed for '%s'.\nRequested size: %ld byte(s) [%.2lf MB]. Total allocated : %ld byte(s) [%.2lf MB]\n", msg, size, (double)size/(double)U_MB, __u_tot_size, (double)__u_tot_size/(double)U_MB);
    }
__u_tot_size += size;

return ptr;
}
int * imalloc(long size, const char * msg)
{
int * ptr;

ptr = (int * )umalloc(sizeof(int)*size, msg);
return ptr;
}


