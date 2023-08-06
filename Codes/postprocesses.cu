#include "kernels.cuh"
#include <cub/cub.cuh>
//#include "ear.h"
#define MAX_MEMORY_SIZE 10000000000
#define MIN(a,b) (a<b?a:b)
typedef unsigned long long ull;
			
void postprocess(int g_n, int g_m, int*& R_d, int*& C_d, int*& F_d, float*& bc_d,int*& d_ld,int* &d_rd, int* &sigma_ld, int*& sigma_rd, 
int *&S_ld,int *&S_rd, int lmemsize,ear_premble **&ear_info, int number_of_SMs,int max_threads_per_block,size_t pitch_ld, 
size_t pitch_rd,size_t pitch_lsigma,size_t pitch_rsigma,size_t pitch_lS,size_t pitch_rS,int start_pos_offset, 
std::vector<int>&level_order,int rmemsize, std::vector<int>&position, double &postproc_time, int *&ear_R, int *&ear_C, 
thrust::device_vector<int> &reachval_d, double &g_active_time)
{
	std::vector<int> local_left_dist,local_right_dist,local_nodes,left_nodes,right_nodes;
	//float *tmp_bc_gpu = new float[g_n];
	//std::cout<<" Ear_Info "<<std::endl;
	for(int i=start_pos_offset;i<(start_pos_offset + lmemsize);i++)
	{
		int v = level_order[i];
		for(int j=ear_R[v];j<ear_R[v+1];j++)
		{
			if(ear_info[j]!=NULL)
			{
				int cnt=0;
				while(cnt<ear_info[j]->cnt)
				{
					if(ear_info[j]->ld[cnt]==0)
						ear_info[j]->ld[cnt] = 1;
					if(ear_info[j]->rd[cnt]==0)
						ear_info[j]->rd[cnt] = 1;

					local_left_dist.push_back(ear_info[j]->ld[cnt]);
					local_right_dist.push_back(ear_info[j]->rd[cnt]);
					local_nodes.push_back(ear_info[j]->nodes[cnt]);
					left_nodes.push_back(v);
					right_nodes.push_back(ear_C[j]);
					cnt++;
				}	
			}
		}
	}
	//std::cout<<" size of local_nodes ====== "<<local_nodes.size()<<" nodes are "<<g_n<<std::endl;	

	if(local_nodes.size()==0)
		return;

	int t_memsize = rmemsize + lmemsize;
	int t_vertices = g_n;
	unsigned long long  int mem_used = ((unsigned long long)rmemsize + lmemsize)*((unsigned long long)g_n)*(4*5);
	unsigned long long int mem_remain = MAX_MEMORY_SIZE - mem_used;
	int allocatable_chunk = mem_remain/((number_of_SMs)*4*5);
	//float *bc_h = new float[g_n];

	int *key_buff_out =  new int[g_n * number_of_SMs]; 
 	int *value_buff_out = new int[g_n * number_of_SMs]; 


	if(allocatable_chunk < number_of_SMs)
	{
		std::cout<<"Not Enough Memory availabel "<<std::endl;
		exit(-1);		
	}	

	//printf("rmemsize = %d lmemsize = %d mem_used = %u \n",rmemsize,lmemsize,mem_used );
	//printf("mem_remain is %llu\n",mem_remain );

	bool *visited_d;
	int *sigma_d, *res_sigma_d;
	int *d_d, *S_d, *key_dist_d, *value_nodes_d, *res_dist_d, *dummy_res_dist_d;
	int *res_seq_d, *dummy_res_dist_out, *res_seq_out, *offset_Q;
	float *delta_d;
	size_t pitch_d, pitch_sigma, pitch_S, pitch_visited, pitch_delta, pitch_key_dist_d, pitch_value_nodes_d, pitch_res_seq_out;
	size_t pitch_dummy_res_dist, pitch_res_seq, pitch_offset_Q, pitch_res_dist, pitch_res_sigma, pitch_dummy_res_dist_out;
	float *tmp_bc_d;
	dim3 dimBlock, dimGrid;
	double total_time=0;	
	GpuTimer timer;	
	//max_threads_per_block = 4;
	allocatable_chunk = (local_nodes.size() < allocatable_chunk ) ? local_nodes.size():allocatable_chunk;
	//printf("allocatable_chunk is %d\n",allocatable_chunk);
	
	dimGrid.x = number_of_SMs;
	dimGrid.y = 1;
	dimGrid.z = 1;

	dimBlock.x = max_threads_per_block;
	dimBlock.y = 1;
	dimBlock.z = 1;

	thrust::device_vector<int>local_left_dist_d(local_nodes.size());
	thrust::device_vector<int>local_right_dist_d(local_nodes.size());
	thrust::device_vector<int>local_nodes_d(local_nodes.size());
	thrust::device_vector<int>left_nodes_d(local_nodes.size());
	thrust::device_vector<int>right_nodes_d(local_nodes.size());
	thrust::device_vector<int>position_d(g_n);

	int n = g_n;
	double pG_time;
	cudaEvent_t psstart,peend;
	thrust::copy(position.begin(), position.end(), position_d.begin());
	
	checkCudaErrors(cudaMallocPitch((void**)&delta_d, &pitch_delta,sizeof(float)*g_n, number_of_SMs));
	checkCudaErrors(cudaMallocPitch(&offset_Q, &pitch_offset_Q,sizeof(int)*g_n, number_of_SMs));
	checkCudaErrors(cudaMallocPitch((void**)&res_dist_d, &pitch_res_dist,sizeof(int)*g_n, number_of_SMs));
	checkCudaErrors(cudaMallocPitch((void**)&res_sigma_d, &pitch_res_sigma,sizeof(int)*g_n, number_of_SMs));
	
	checkCudaErrors(cudaMalloc(&dummy_res_dist_d,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMalloc(&res_seq_d,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMalloc(&dummy_res_dist_out,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMalloc(&res_seq_out,sizeof(int)*g_n*number_of_SMs));

	//checkCudaErrors(cudaMalloc(&tmp_bc_d,sizeof(float)*g_n));

	checkCudaErrors(cudaMemset(dummy_res_dist_out,0,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMemset(res_seq_out,0,sizeof(int)*g_n*number_of_SMs));

	local_left_dist_d = local_left_dist;
	local_right_dist_d = local_right_dist;
	local_nodes_d = local_nodes;
	left_nodes_d = left_nodes;
	right_nodes_d = right_nodes;

	int num_iterations = (int)ceil(((double)local_nodes.size())/number_of_SMs);

	for(int i=0;i < num_iterations;i++)
	{
		int limit = MIN(number_of_SMs+(i*number_of_SMs),local_nodes.size());//local_nodes.size();
		int begin = (i*number_of_SMs);
		//checkCudaErrors(cudaMemset(tmp_bc_d,0,sizeof(float)*g_n));
		timer.Start();			

		bc_postprocess1<<<dimGrid,dimBlock>>>(bc_d, sigma_ld, sigma_rd, d_ld, d_rd, pitch_lsigma, pitch_rsigma, pitch_ld, pitch_rd, 
		thrust::raw_pointer_cast(local_left_dist_d.data()),	thrust::raw_pointer_cast(local_right_dist_d.data()), 
		thrust::raw_pointer_cast(local_nodes_d.data()), thrust::raw_pointer_cast(left_nodes_d.data()), thrust::raw_pointer_cast(right_nodes_d.data()),
		R_d,C_d,limit, g_n, thrust::raw_pointer_cast(position_d.data()), lmemsize, start_pos_offset, res_dist_d, pitch_res_dist, 
		res_sigma_d, pitch_res_sigma, dummy_res_dist_d, res_seq_d, offset_Q, pitch_offset_Q, g_m, F_d, 
		thrust::raw_pointer_cast(reachval_d.data()),begin);

		checkCudaErrors(cudaPeekAtLastError());
		timer.Stop();
		total_time += timer.Elapsed();	

		for(int j=0;j<number_of_SMs;j++)
		{
			void  *d_temp_storage = NULL;
			size_t temp_storage_bytes = 0;
			cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, dummy_res_dist_d + (j * g_n), dummy_res_dist_out + (j * g_n), res_seq_d + (j * g_n), res_seq_out + (j * g_n), g_n);
			checkCudaErrors(cudaDeviceSynchronize());
			//d_temp_storage = new  int[temp_storage_bytes];
			cudaMalloc(&d_temp_storage, temp_storage_bytes);
			cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, dummy_res_dist_d + (j * g_n), dummy_res_dist_out + (j * g_n), res_seq_d + (j * g_n), res_seq_out + (j * g_n), g_n);
			checkCudaErrors(cudaDeviceSynchronize());
			//delete[] d_temp_storage;	
			checkCudaErrors(cudaFree(d_temp_storage));
		}
		
		timer.Start();			
		bc_postprocess2<<<dimGrid,dimBlock>>>(bc_d, R_d, C_d, g_n, delta_d, pitch_delta, res_dist_d, res_sigma_d, res_seq_out,
		dummy_res_dist_out, pitch_res_dist, pitch_res_sigma, offset_Q, pitch_offset_Q, g_m, F_d, 
		thrust::raw_pointer_cast(reachval_d.data()), limit,begin, thrust::raw_pointer_cast(local_nodes_d.data()));

		checkCudaErrors(cudaPeekAtLastError());

		timer.Stop();
		total_time += timer.Elapsed();	

	
		//checkCudaErrors(cudaMemcpy(bc_h,bc_d,sizeof(float)*g_n,cudaMemcpyDeviceToHost));
		/*checkCudaErrors(cudaMemcpy(tmp_bc_gpu,tmp_bc_d,sizeof(float)*g_n,cudaMemcpyDeviceToHost));
		printf("TEMP PArtial BC after Postprocess Outside for loop\n");	
		for(int h=0; h<g_n; h++)
		{
			printf("%d = %f\n",h, tmp_bc_gpu[h] );
		}*/

	}

	local_nodes_d.clear();
	local_nodes_d.shrink_to_fit();	
	left_nodes_d.clear();
	left_nodes_d.shrink_to_fit();
	right_nodes_d.clear();
	right_nodes_d.shrink_to_fit();
	local_left_dist_d.clear();
	local_left_dist_d.shrink_to_fit();
	local_right_dist_d.clear();
	local_right_dist_d.shrink_to_fit();	
	position_d.clear();
	position_d.shrink_to_fit();
	
	checkCudaErrors(cudaFree(delta_d));
	checkCudaErrors(cudaFree(offset_Q));
	checkCudaErrors(cudaFree(res_dist_d));
	checkCudaErrors(cudaFree(dummy_res_dist_d));
	checkCudaErrors(cudaFree(res_sigma_d));
	checkCudaErrors(cudaFree(res_seq_d));
	checkCudaErrors(cudaFree(dummy_res_dist_out));
	checkCudaErrors(cudaFree(res_seq_out));
	
	printf("postprocess time is %f\n",total_time);
	g_active_time += total_time;

}

__device__ void check_neighbours(const int *&R, const int *&C, int *&d_row, int *&dummy_d_row, int *&sigma_row, int &node, 
int &left, int &right)
{
		int prv_node = node;
		int a = C[R[node]];
		int b = C[R[node] + 1];
		int dist = 1;
			
		//while(a!=left && a!=right)
		while((R[a+1]-R[a])==2)
		{
		//	printf("aaaa= %d node %d left %d right %d\n",a,node, left,right);
			if(d_row[a] > dist)
			{
				d_row[a] = dist;
				dummy_d_row[a] = d_row[a];
				sigma_row[a] = 1;
			}
			else 
			{
				if(d_row[a] == dist)
					sigma_row[a] += 1;	
			}
			
			++dist;

			if(C[R[a]]!=prv_node)
			{
				prv_node = a;
				a = C[R[a]];
			}else
			{
				prv_node = a;
				a = C[R[a] + 1];
			}	
			//printf("aaa= %d node %d left %d right %d\n",a,node, left,right);
		}	

		dist = 1;
		prv_node = node;
		//while(b!=right && b!=left)
		while((R[b+1]-R[b])==2)
		{
		//	printf("bbbb= %d node %d left %d right %d\n",b,node, left,right);
			if(d_row[b] > dist)
			{
				d_row[b] = dist;
				dummy_d_row[b] = d_row[b];
				sigma_row[b] = 1;
			}
			else 
			{
				if(d_row[b] == dist)
					sigma_row[b] += 1;	
			}
			
			++dist;

			if(C[R[b]]!=prv_node)
			{
				prv_node = b;
				b = C[R[b]];
			}else
			{
				prv_node = b;	
				b = C[R[b] + 1];
			}	
			//printf("bbb= %d node %d left %d right %d\n",b,node, left,right);
		}
}

__global__ void bc_postprocess1(float* bc, int *sigma_ld,int *sigma_rd, int* d_ld, int* d_rd, size_t pitch_lsigma, size_t pitch_rsigma, 
	size_t pitch_ld,size_t pitch_rd, int *left_dist, int *right_dist, int *local_nodes, int *left_nodes, int *right_nodes, 
	const int *R_d, const int *C_d,int limit, const int n, int *position,int lmemsize, int start_pos_offset, int *res_dist,
	size_t pitch_res_dist, int *res_sigma, size_t pitch_res_sigma, int *dummy_res_dist, int *res_seq, 
	int *offset_Q, size_t pitch_offset_Q, const int m,  const int *F, int *reachval, int begin)
{
	int tid = threadIdx.x;
	int index = 0;
	int *d_lrow, *d_rrow;//*S_lrow, *S_rrow;
	int *sigma_lrow, *sigma_rrow;

	int *dummy_res_dist_row = (int*)(dummy_res_dist+ blockIdx.x * n);
	
	int *res_seq_row = (int*)(res_seq + blockIdx.x * n);
	int *res_dist_row = (int*)((char*)res_dist + blockIdx.x * pitch_res_dist);
	int *res_sigma_row = (int*)((char*)res_sigma + blockIdx.x * pitch_res_sigma);	

	__syncthreads();

	for(int z=begin+blockIdx.x; z<limit; z+=gridDim.x)
	{	
		//if((z * gridDim.x + blockIdx.x) < limit)
		{	
			for(int k=threadIdx.x; k < n; k+=blockDim.x)
			{
				dummy_res_dist_row[k] = INT_MAX;
				//delta_row[k] = 0;	
				res_dist_row[k] = INT_MAX;
			}

			int ldz = left_dist[z];
			int rdz = right_dist[z];
			int current_node = local_nodes[z];
			int left = left_nodes[z];
			int right = right_nodes[z];

			d_lrow = (int*)((char*)d_ld + (position[left] - start_pos_offset) * pitch_ld);
			sigma_lrow = (int*)((char*)sigma_ld + (position[left] - start_pos_offset) * pitch_lsigma);
			//S_lrow = (int*)((char*)S_ld + (position[left] - start_pos_offset) * pitch_lS);
			
			if(position[right] < (start_pos_offset + lmemsize))
			{
				d_rrow = (int*)((char*)d_ld + (position[right] - start_pos_offset) * pitch_ld);
				sigma_rrow = (int*)((char*)sigma_ld + (position[right] - start_pos_offset) * pitch_lsigma);
				//S_rrow = (int*)((char*)S_ld + (position[right] - start_pos_offset) * pitch_lS);	
			}else
			{
				d_rrow = (int*)((char*)d_rd + (position[right] - (start_pos_offset + lmemsize) ) * pitch_rd);
				sigma_rrow = (int*)((char*)sigma_rd + (position[right] - (start_pos_offset + lmemsize) ) * pitch_rsigma);
				//S_rrow = (int*)((char*)S_rd + (position[right] - (start_pos_offset + lmemsize) ) * pitch_rS);		
			}
			__syncthreads();

			for(int k=threadIdx.x; k<n; k+=blockDim.x)
			{
				int dlk = d_lrow[k];
				int drk = d_rrow[k];

				if((dlk+ ldz) < (drk + rdz ))																												
				{
					res_dist_row[k] = dlk + ldz;
					res_sigma_row[k] = sigma_lrow[k];
				}
				else if((dlk+ ldz) > (drk + rdz ))	
				{
					res_dist_row[k] = drk + rdz;
					res_sigma_row[k] = sigma_rrow[k];
				}
				else
				{
					res_dist_row[k] = drk + rdz;
					res_sigma_row[k] = sigma_rrow[k] + sigma_lrow[k];
				}
				dummy_res_dist_row[k] = res_dist_row[k];	
				res_seq_row[k] = k;
			}	
			__syncthreads();
			
			if(threadIdx.x==0)
			{
				dummy_res_dist_row[current_node] = 0;
				res_dist_row[current_node] = 0;	
				res_sigma_row[current_node] = 1;
				check_neighbours(R_d, C_d, res_dist_row, dummy_res_dist_row, res_sigma_row, current_node, left, right);
			}
		}	
	}
}


__global__ void bc_postprocess2(float *bc, const int *R_d, const int *C_d, int n, float *delta_d, size_t pitch_delta, int *res_dist, 
int *res_sigma, int *res_seq_out, int *dummy_res_dist_out, size_t pitch_res_dist, size_t pitch_res_sigma,
int *offset_Q, size_t pitch_offset_Q, int m, const int *F, int *reachval, int limit,int begin, int *local_nodes)
{

	int tid = threadIdx.x;

	for(int z = begin+blockIdx.x; z < limit; z+=gridDim.x)
	{
		int current_node = local_nodes[z];
		int *dummy_res_dist_out_row = (int *)(dummy_res_dist_out + blockIdx.x * n);
		int *S_row = (int*)(res_seq_out + blockIdx.x*n);
		
		int *res_dist_row = (int*)((char*)res_dist + blockIdx.x*pitch_res_dist);
		int *res_sigma_row = (int*)((char*)res_sigma + blockIdx.x*pitch_res_sigma);
		float *delta_row = (float*)((char*)delta_d + blockIdx.x * pitch_delta);
		
		__shared__ int *offset_Q_row;
		__shared__ int offset_ends;
		__shared__ int current_depth;
		__shared__ int start;
		
		__syncthreads();

		if(tid==0)
		{
			offset_Q_row = (int *)((char*)offset_Q + pitch_offset_Q * blockIdx.x);
			offset_Q_row[0] = 0;
			offset_Q_row[ dummy_res_dist_out_row[n-1] + 1] = n;  
			offset_ends = dummy_res_dist_out_row[n-1] + 1;
			current_depth = offset_ends - 1;
			start = S_row[0];
		}
		
		__syncthreads();
		
		for(int i=threadIdx.x; i<n; i+=blockDim.x)
		{
			if(i!=0)
			{
				if(dummy_res_dist_out_row[i]!=dummy_res_dist_out_row[i-1])
				{
					offset_Q_row[ dummy_res_dist_out_row[i] ] = i; 
				}
			}
			delta_row[i] = 0;	
		}

		__syncthreads();
		
		while(current_depth > 0)
		{
			int stack_iter_len = offset_Q_row[current_depth+1]-offset_Q_row[current_depth];
			if(stack_iter_len>512)
			{
				for(int kk=threadIdx.x; kk<2*m; kk+=blockDim.x)
				{
					int w = F[kk];
					if(res_dist_row[w] == current_depth)
					{
						int v = C_d[kk];
						if(res_dist_row[v] == (res_dist_row[w]+1))
						{
							float change = (res_sigma_row[w]/(float)res_sigma_row[v])*(1.0f * (reachval[v]+1) + delta_row[v]);
							atomicAdd(&delta_row[w],change);
						}		
					}
				}
			}
			else 
			{
				for(int kk=threadIdx.x+offset_Q_row[current_depth]; kk<offset_Q_row[current_depth+1]; kk+=blockDim.x)
				{
					int w = S_row[kk];
					float dsw = 0;
					float sw = (float)res_sigma_row[w];
					for(int bb=R_d[w]; bb<R_d[w+1]; bb++)
					{
						int v = C_d[bb];
						if(res_dist_row[v] == (res_dist_row[w]+1))
						{
							dsw += (sw/(float)res_sigma_row[v])*(1.0f * (reachval[v]+1) + delta_row[v]);
							//if(current_node==2)
							//		printf("current node = %d parent %d child %d reachval %d start %d current_depth %d sigma[%d] %u dsw %f t_nodes %d\n",current_node, w,v,reachval[v],current_node,current_depth,v,res_sigma_row[v],dsw,n);
						}
					}
					delta_row[w] += dsw ;
				}
			}
			__syncthreads();
			if(tid == 0)
			{
				current_depth--;
			}
			__syncthreads();
		}
		for(int kk=threadIdx.x; kk<n; kk+=blockDim.x)
		{
			atomicAdd(&bc[kk], (delta_row[kk] * (reachval[current_node]+1))); //Would need to check that kk != i here, but delta_row[kk] is guaranteed to be 0.
			//atomicAdd(&tmp_bc[kk],(delta_row[kk] * (reachval[current_node]+1)));
		}
		//==========================
	}

}


	