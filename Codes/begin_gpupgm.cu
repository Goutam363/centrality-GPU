#include <set>
#include <vector>
#include "kernels.cuh"
#define THREADS_PER_BLOCK 1024
#define WARP_SIZE 32
using namespace std;

void choose_device(int &max_threads_per_block, int &number_of_SMs)
{
	int count;
	checkCudaErrors(cudaGetDeviceCount(&count));
	cudaDeviceProp prop;

	//if(op.device == -1)
	{
		int maxcc=0, bestdev=0;
		for(int i=0; i<count; i++)
		{
			checkCudaErrors(cudaGetDeviceProperties(&prop,i));
			if((prop.major + 0.1*prop.minor) > maxcc)
			{
				maxcc = prop.major + 0.1*prop.minor;
				bestdev = i;
			}	
		}

		checkCudaErrors(cudaSetDevice(bestdev));
		checkCudaErrors(cudaGetDeviceProperties(&prop,bestdev));
	}
	
	//std::cout << "Chosen Device: " << prop.name << std::endl;
	//std::cout << "Compute Capability: " << prop.major << "." << prop.minor << std::endl;
	//std::cout << "Number of Streaming Multiprocessors: " << prop.multiProcessorCount << std::endl;
	//std::cout << "Size of Global Memory: " << prop.totalGlobalMem/(float)(1024*1024*1024) << " GB" << std::endl << std::endl;
	//std::cout << "MaxThreadPerBolck: " << prop.maxThreadsPerBlock<< std::endl << std::endl;
	//std::cout << "MultiProcessorCount: " << prop.multiProcessorCount<< std::endl << std::endl;
	max_threads_per_block = prop.maxThreadsPerBlock;
	number_of_SMs = prop.multiProcessorCount;
}


void begin_gpu(int *&R, int *&C, int *&F, int n, int m, ear_premble **&ear_info, ear_premble**&self_ear, std::vector<int> &req_nodes, 
std::vector<int> &artnodes, float *&bc_gpu, double &g_free_time)
{

	int max_threads_per_block, number_of_SMs;
	std::set<int> source_vertices;

	choose_device(max_threads_per_block,number_of_SMs);

	//if( n < THREADS_PER_BLOCK && n > WARP_SIZE)
	if( n < THREADS_PER_BLOCK )
	{
		max_threads_per_block = n;
	}
	/*else
	{
		if(n < WARP_SIZE)
		{
			return ;
		}	
	}*/	
	
	if(req_nodes.size() < number_of_SMs)
		number_of_SMs = req_nodes.size();

	bc_gpu_free(n, m/2, R, C, F, max_threads_per_block, number_of_SMs, source_vertices, req_nodes, artnodes, bc_gpu, g_free_time);
	
	std::cout<<"Free vertices processing finished"<<std::endl;

}

void begin_gpu(int *&R, int *&C, int *&F, int n, int m, ear_premble **&ear_info, std::vector<int>&parents, std::vector<int>&levels, std::vector<int>&level_order, 
std::vector<int>&level_offset, int max_earlevel, std::vector<int>&artnodes, float *&bc_gpu, double&g_free_time, int *&ear_R, int *&ear_C)
{

	int max_threads_per_block, number_of_SMs;
	std::set<int> source_vertices;

	choose_device(max_threads_per_block,number_of_SMs);

	if( n < THREADS_PER_BLOCK )
	{
		max_threads_per_block = n;
	}

	bc_gpu_active(n, m/2, R, C, F, max_threads_per_block, number_of_SMs, source_vertices, parents, levels, level_order, level_offset, 
	ear_info, max_earlevel, bc_gpu, g_free_time, ear_R, ear_C, artnodes);
}