#ifndef GREEDY_MATCHER
#define GREEDY_MATCHER
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "constants.cuh"

#include "matchgpu.h"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/zip_function.h>
#include <thrust/unique.h>
struct is_less_than_0
{
  __host__ __device__ int operator()(int &x)
  {
    if (x > -1)
      return x;
    else
      return -1;
  }
};

struct GreedyMatcher
{

  GreedyMatcher(int _n, int *_cmatch, int *_rmatch):
  n(_n),cmatch(_cmatch),rmatch(_rmatch)
  {
    m_h.resize(n);
    c_h.resize(1, 1);
    m_d.resize(n, 0);
    req_d.resize(n, 0);
    c_d.resize(1, 1);
    cudaMallocHost((void**)&c_Pinned, sizeof(int)); // host pinned
  }
  int maxMatch(int * match)
  {
    //int *m_d_ptr = thrust::raw_pointer_cast(m_d.data());
    int *m_d_ptr = thrust::raw_pointer_cast(m_d.data());
    int *req_d_ptr = thrust::raw_pointer_cast(req_d.data());
    int *c_ptr = thrust::raw_pointer_cast(c_d.data());
    srand(1);
    /*computing MM */
    int matchround = 0;
    int dimGrid = (n + THREADS_PER_BLOCK) / THREADS_PER_BLOCK;
    memset(c_Pinned, 1, sizeof(int));
    while (c_Pinned[0] && ++matchround < NR_MAX_MATCH_ROUNDS)
    {
      //printf("match round %d n %d\n", matchround,n);
      gaSelect<<<dimGrid, THREADS_PER_BLOCK>>>(m_d_ptr, c_ptr, n, rand());
      cudaMemcpy(c_Pinned, c_ptr, sizeof(int), cudaMemcpyDeviceToHost);
      
      grRequestEdgeList<<<dimGrid, THREADS_PER_BLOCK>>>(rmatch, cmatch, req_d_ptr, m_d_ptr, c_ptr,n);
      
      grRespondEdgeList<<<dimGrid, THREADS_PER_BLOCK>>>(rmatch, cmatch, req_d_ptr, m_d_ptr,n);
      
      gMatchEdgeList<<<dimGrid, THREADS_PER_BLOCK>>>(rmatch, cmatch, m_d_ptr,req_d_ptr, n);
      gTentativelyKill<<<dimGrid, THREADS_PER_BLOCK>>>(m_d_ptr,req_d_ptr,n);
      gRestoreLife<<<dimGrid, THREADS_PER_BLOCK>>>(rmatch, cmatch, m_d_ptr,req_d_ptr, n);
    }

    using namespace thrust::placeholders;
    thrust::for_each(m_d.begin(), m_d.end(), _1 -= 4);
    thrust::transform(m_d.begin(), m_d.end(), m_d.begin(), is_less_than_0()); // in-place transformation
    //m_h = m_d;
    //for (int i = 0; i < m_h.size(); ++i)
    //    printf("%d %d\n", i, m_h[i]);
    //for (int i = 0; i < m_h.size(); ++i)
    //  mate[i] = m_h[i];
    cudaMemcpy(match, m_d_ptr, sizeof(int)*(m_d.size()), cudaMemcpyDeviceToHost);

    int numAugmented = thrust::count_if(m_d.begin(), m_d.end(), _1 > -1);
    return numAugmented / 2;
  }

  int n;
  int * cmatch, *rmatch;
  thrust::host_vector<int> m_h;
  thrust::device_vector<int> c_h;

  thrust::device_vector<int> m_d;
  thrust::device_vector<int> req_d;
  thrust::device_vector<int> c_d;
  int * c_Pinned;
};
#endif