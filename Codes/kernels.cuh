//#ifndef BC_KERNELS
//#define BC_KERNELS
#ifndef BC_UTIL
#define BC_UTIL

#include <vector>
#include <cuda.h>
#include <set>
#include <cstdlib>
#include <iostream>
//#include <cuda.h>
//#include <getopt.h>
//#include <cmath>

#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include "ear_graph.h"
#include "utils.h"


#ifndef checkCudaErrors
#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

// These are the inline versions for all of the SDK helper functions
inline void __checkCudaErrors(cudaError_t err, const char *file, const int line)
{   
    if(cudaSuccess != err)
    {   
        std::cerr << "CUDA Error = " << err << ": " << cudaGetErrorString(err) << " from file " << file  << ", line " << line << std::endl;
        exit(EXIT_FAILURE);
    }
}
#endif


__device__ void bitonic_sort_free(int *values, int N);
//__device__ void bitonic_sort_active(int *values, int N);
//__device__ void bitonic_sort(int *values, int N);



void bc_gpu_free(int g_nodes, int g_edges, int *&g_R, int *&g_C, int *&g_F, int max_threads_per_block, int number_of_SMs,
const std::set<int> &source_verticess, std::vector<int> &req_nodes_hs, std::vector<int> &reachval_hs, float *&bc_gpus, double &g_free_time);

void bc_gpu_active(int n, int m, int *&R, int *&C, int *&F, int max_threads_per_block, int number_of_SMs, const std::set<int> &source_vertices,
std::vector<int>&parents, std::vector<int>&levels, std::vector<int>&level_order, std::vector<int>&level_offset, ear_premble **&ear_info, 
int max_earlevel, float*&bc_gpu, double&g_active_time, int *&ear_R, int *&ear_C, std::vector<int>&reachval_h);

void postprocess(int g_n, int g_m, int*& R_d, int*& C_d, int*& F_d, float*& bc_d,int*& d_ld,int* &d_rd, int* &sigma_ld, int*& sigma_rd, 
int *&S_ld,int *&S_rd, int lmemsize,ear_premble **&ear_info, int number_of_SMs,int max_threads_per_block,size_t pitch_ld, 
size_t pitch_rd,size_t pitch_lsigma,size_t pitch_rsigma,size_t pitch_lS,size_t pitch_rS,int start_pos_offset, 
std::vector<int>&level_order,int rmemsize, std::vector<int>&position, double &postproc_time, int *&ear_R, int *&ear_C, 
thrust::device_vector<int> &reachval_d, double &g_active_time);

__global__ void bc_kgpu_free(float *bc, const int *R, const int *C, const int *F, const int n, const int m,int*d, 
unsigned long long *sigma,float *delta, int *Q, int *Q2,int*S,int *endpoints, int *next_source,size_t pitch_d, size_t pitch_sigma,
size_t pitch_delta, size_t pitch_Q, size_t pitch_Q2,size_t pitch_S, size_t pitch_endpoints, int start, int end, int *jia, int *diameters, 
int *source_vertices, bool approx, int* req_nodes, int *artnodes);

__global__ void bc_kgpu_active(float *bc, const int *R, const int *C, const int *F, const int n, const int m,int*d, int *ld,
int *rd,int *sigma,int *lsigma,int *rsigma, float *delta, int *Q, int *Q2,int*S,
int *lS,int *rS,int *endpoints, int *next_source,size_t pitch_d, size_t pitch_ld,size_t pitch_rd,size_t pitch_sigma,
size_t pitch_lsigma,size_t pitch_rsigma,size_t pitch_delta, size_t pitch_Q, size_t pitch_Q2,size_t pitch_S,size_t pitch_lS,
size_t pitch_rS,size_t pitch_endpoints, int start, int end, int *jia, int *diameters,
int *active_vertices, int rmemsize, int *reachval);

__global__ void bc_postprocess1(float* bc, int *sigma_ld,int *sigma_rd, int* d_ld, int* d_rd, size_t pitch_lsigma, size_t pitch_rsigma, 
size_t pitch_ld,size_t pitch_rd, int *left_dist, int *right_dist, int *local_nodes, int *left_nodes, int *right_nodes, 
const int *R_d, const int *C_d,int limit, const int n, int *position,int lmemsize, int start_pos_offset, int *res_dist,
size_t pitch_res_dist, int *res_sigma, size_t pitch_res_sigma, int *dummy_res_dist, int *res_seq, 
int *offset_Q, size_t pitch_offset_Q, const int m,  const int *F, int *reachval, int begin);

__global__ void bc_postprocess2(float *bc, const int *R_d, const int *C_d, int n, float *delta_d, size_t pitch_delta, int *res_dist, 
int *res_sigma, int *res_seq_out, int *dummy_res_dist_out, size_t pitch_res_dist, size_t pitch_res_sigma,
int *offset_Q, size_t pitch_offset_Q, int m, const int *F, int *reachval, int limit,int begin, int *local_nodes);

__device__ void check_neighbours(int *&R, int *&C, int *&d_row, int *&dummy_d_row, int *&sigma_row, int &node,
 int &left, int &right);

#endif