#include "kernels.cuh"
#include <cub/cub.cuh>
//#include "ear.h"
#define MAX_MEMORY_SIZE 10000000000
#define MIN(a,b) (a<b?a:b)
typedef unsigned long long ull;

//int lnode,int rnode,int cnt,int td,int *ld,int *rd,int *nodes

//void postprocess(float* bc_d,int* d_ld,int *d_rd,int* sigma_ld,int* sigma_rd,int *S_ld,int *S_rd,graph gnew,std::vector<int> reordering,int lmemsize,std::vector<bool>required_storage,dg2** ear_info,std::vector<int>reordering_map,thrust::device_vector<int>&nomap_d,int number_of_SMs,int max_threads_per_block,size_t pitch_ld,size_t pitch_rd,size_t pitch_lsigma,size_t pitch_rsigma,size_t pitch_lS,size_t pitch_rS)
void postprocess(int g_n, int g_m, int*& R_d, int*& C_d, int*& F_d, float*& bc_d,int*& d_ld,int* &d_rd, unsigned long long* &sigma_ld, unsigned long long*& sigma_rd, 
int *&S_ld,int *&S_rd, int lmemsize,ear_premble **&ear_info, int number_of_SMs,int max_threads_per_block,size_t pitch_ld, 
size_t pitch_rd,size_t pitch_lsigma,size_t pitch_rsigma,size_t pitch_lS,size_t pitch_rS,int start_pos_offset, 
std::vector<int>&level_order,int rmemsize, std::vector<int>&position, double &postproc_time, int *&ear_R, int *&ear_C, 
thrust::device_vector<int> &reachval_d, double &g_active_time)
{
	std::vector<int> local_left_dist,local_right_dist,local_nodes,left_nodes,right_nodes;
	float *tmp_bc_gpu = new float[g_n];
	std::cout<<" Ear_Info "<<std::endl;
	for(int i=start_pos_offset;i<(start_pos_offset + lmemsize);i++)
	{
		int v = level_order[i];
		for(int j=ear_R[v];j<ear_R[v+1];j++)
		{
			if(ear_info[j]!=NULL)
			{
				int cnt=0;
				//if(ear_info[j]->cnt > 1)
					std::cout<<v<<" "<<ear_C[j]<<std::endl;
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
				//	if(ear_info[j]->cnt > 1)
						std::cout<<ear_info[j]->ld[cnt]<<" "<<ear_info[j]->nodes[cnt]<<" "<<ear_info[j]->rd[cnt]<<std::endl;
					cnt++;
				}	
			}
		}
	}
	std::cout<<" size of local_nodes ====== "<<local_nodes.size()<<" nodes are "<<g_n<<std::endl;	

	if(local_nodes.size()==0)
		return;

	int t_memsize = rmemsize + lmemsize;
	int t_vertices = g_n;
	unsigned long long  int mem_used = ((unsigned long long)rmemsize + lmemsize)*((unsigned long long)g_n)*(4*5);
	unsigned long long int mem_remain = MAX_MEMORY_SIZE - mem_used;
	int allocatable_chunk = mem_remain/((number_of_SMs)*4*5);
	float *bc_h = new float[g_n];

	int *key_buff_out =  new int[g_n * number_of_SMs]; 
 	int *value_buff_out = new int[g_n * number_of_SMs]; 


	if(allocatable_chunk < number_of_SMs)
	{
		std::cout<<"Not Enough Memory availabel "<<std::endl;
		exit(-1);		
	}	

//	std::cout<<"mem_used "<<mem_used<<" "<<((rmemsize + lmemsize)*g_n*4*5)<<" dummy_mem_used "<<dummy_mem_used<<" int_dummy_mem_used "<<int_dummy_mem_used<<std::endl;
//	std::cout<<"mem_used "<<mem_used<<" mem_remain "<<mem_remain<<" allocatable_chunk "<<allocatable_chunk<<std::endl;
	printf("rmemsize = %d lmemsize = %d mem_used = %u \n",rmemsize,lmemsize,mem_used );
	printf("mem_remain is %llu\n",mem_remain );
//	printf("allocatable_chunk is %d\n",allocatable_chunk);

	//bool *visited_d;
	unsigned long long *sigma_d, *res_sigma_d;
	int *d_d, *S_d, *key_dist_d, *value_nodes_d, *res_dist_d, *dummy_res_dist_d;
	int *res_seq_d, *dummy_res_dist_out, *res_seq_out, *offset_Q;
	float *delta_d;
	size_t pitch_d, pitch_sigma, pitch_S, pitch_delta, pitch_key_dist_d, pitch_value_nodes_d, pitch_res_seq_out;
	size_t pitch_dummy_res_dist, pitch_res_seq, pitch_offset_Q, pitch_res_dist, pitch_res_sigma, pitch_dummy_res_dist_out, pitch_visited;
	float *tmp_bc_d;
	bool *visited_d;
	dim3 dimBlock, dimGrid;
	double total_time=0;	
	GpuTimer timer;	
	//max_threads_per_block = 4;
	allocatable_chunk = (local_nodes.size() < allocatable_chunk ) ? local_nodes.size():allocatable_chunk;
	printf("allocatable_chunk is %d\n",allocatable_chunk);
	
	dimGrid.x = number_of_SMs;
	dimGrid.y = 1;
	dimGrid.z = 1;

	dimBlock.x = max_threads_per_block;
	dimBlock.y = 1;
	dimBlock.z = 1;

	thrust::device_vector<int>local_left_dist_d(number_of_SMs);
	thrust::device_vector<int>local_right_dist_d(number_of_SMs);
	thrust::device_vector<int>local_nodes_d(number_of_SMs);
	thrust::device_vector<int>left_nodes_d(number_of_SMs);
	thrust::device_vector<int>right_nodes_d(number_of_SMs);
	thrust::device_vector<int>position_d(g_n);

	int n = g_n;
	
	double pG_time;
	cudaEvent_t psstart,peend;
	thrust::copy(position.begin(), position.end(), position_d.begin());
	

	checkCudaErrors(cudaMallocPitch((void**)&delta_d, &pitch_delta,sizeof(float)*g_n, number_of_SMs));
	checkCudaErrors(cudaMallocPitch(&offset_Q, &pitch_offset_Q,sizeof(int)*g_n, number_of_SMs));
	checkCudaErrors(cudaMallocPitch((void**)&res_dist_d, &pitch_res_dist,sizeof(int)*g_n, number_of_SMs));
	checkCudaErrors(cudaMalloc((void**)&dummy_res_dist_d,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMallocPitch((void**)&res_sigma_d, &pitch_res_sigma,sizeof(unsigned long long)*g_n, number_of_SMs));
	checkCudaErrors(cudaMalloc((void**)&res_seq_d,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMalloc(&dummy_res_dist_out,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMalloc(&res_seq_out,sizeof(int)*g_n*number_of_SMs));
	checkCudaErrors(cudaMalloc((void**)&tmp_bc_d,sizeof(float)*g_n));
	checkCudaErrors(cudaMallocPitch((void**)&visited_d,&pitch_visited,sizeof(bool)*g_n,number_of_SMs));

	//for(int i=0;i < local_nodes.size();i+=1)
	for(int i=0;i < local_nodes.size();i+=number_of_SMs)
	{
		int limit = MIN(number_of_SMs,local_nodes.size()-i);
		//int limit = 1;
		thrust::copy(local_left_dist.begin() + i, local_left_dist.begin() + i + limit, local_left_dist_d.begin());
		thrust::copy(local_right_dist.begin() + i, local_right_dist.begin() + i + limit, local_right_dist_d.begin());
		thrust::copy(local_nodes.begin() + i, local_nodes.begin() + i + limit, local_nodes_d.begin());
		thrust::copy(left_nodes.begin() + i, left_nodes.begin() + i + limit, left_nodes_d.begin());
		thrust::copy(right_nodes.begin() + i, right_nodes.begin() + i + limit, right_nodes_d.begin());

		checkCudaErrors(cudaMemset(dummy_res_dist_out,0,sizeof(int)*g_n*number_of_SMs));
		checkCudaErrors(cudaMemset(res_seq_out,0,sizeof(int)*g_n*number_of_SMs));
		checkCudaErrors(cudaMemset(tmp_bc_d,0,sizeof(float)*g_n));
		//printf("Going for bc_postprocess1 limit %d remains %d\n",limit,local_nodes.size()-i);

		timer.Start();			
		bc_postprocess1<<<dimGrid,dimBlock>>>(bc_d, sigma_ld, sigma_rd, d_ld, d_rd, S_ld, S_rd, pitch_lsigma, pitch_rsigma, pitch_ld, pitch_rd, 
		pitch_lS, pitch_rS, thrust::raw_pointer_cast(local_left_dist_d.data()),	thrust::raw_pointer_cast(local_right_dist_d.data()), 
		thrust::raw_pointer_cast(local_nodes_d.data()), thrust::raw_pointer_cast(left_nodes_d.data()), thrust::raw_pointer_cast(right_nodes_d.data()),
		R_d,C_d,limit, g_n, thrust::raw_pointer_cast(position_d.data()), lmemsize, start_pos_offset, delta_d, pitch_delta, res_dist_d,
		pitch_res_dist, res_sigma_d, pitch_res_sigma, dummy_res_dist_d, res_seq_d, visited_d, pitch_visited);

		checkCudaErrors(cudaPeekAtLastError());
		//checkCudaErrors(cudaDeviceSynchronize());
		timer.Stop();
		total_time += timer.Elapsed();	

		for(int k=0; k < number_of_SMs; k++)
		{	
			if(k<limit)
			{
				void  *d_temp_storage = NULL;
				size_t temp_storage_bytes = 0;
				cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, dummy_res_dist_d + (k*g_n), dummy_res_dist_out + (k*g_n), res_seq_d + (k*g_n), res_seq_out + (k*g_n), g_n);
				cudaMalloc(&d_temp_storage, temp_storage_bytes);
				cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, dummy_res_dist_d + (k*g_n), dummy_res_dist_out + (k*g_n), res_seq_d + (k*g_n), res_seq_out + (k*g_n), g_n);
				checkCudaErrors(cudaDeviceSynchronize());
				checkCudaErrors(cudaFree(d_temp_storage));
			}	
		}	
		//printf("Going for bc_postprocess2\n");

		timer.Start();			
		bc_postprocess2<<<dimGrid,dimBlock>>>(bc_d, R_d, C_d, g_n, delta_d, pitch_delta, res_dist_d, res_sigma_d, res_seq_out,
		dummy_res_dist_out, pitch_res_dist, pitch_res_sigma, offset_Q, pitch_offset_Q, g_m, F_d, 
		thrust::raw_pointer_cast(reachval_d.data()), limit, tmp_bc_d);

		checkCudaErrors(cudaPeekAtLastError());

		timer.Stop();
		total_time += timer.Elapsed();	
	
		checkCudaErrors(cudaMemcpy(bc_h,bc_d,sizeof(float)*g_n,cudaMemcpyDeviceToHost));
		//printf("After postprocess\n");

		/*checkCudaErrors(cudaMemcpy(tmp_bc_gpu,tmp_bc_d,sizeof(float)*g_n,cudaMemcpyDeviceToHost));
		printf("PArtial betweenness centrality postprocess \n");
		for(int h=0; h<g_n; h++)
		{
			printf("%d = %f\n",h, tmp_bc_gpu[h]/2 );
		}*/

	}

	printf("Outside for loop\n");
	
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
	
	//printf("Outside for loop free memory 1 \n");	
	checkCudaErrors(cudaFree(delta_d));
	checkCudaErrors(cudaFree(offset_Q));
	checkCudaErrors(cudaFree(res_dist_d));
	checkCudaErrors(cudaFree(dummy_res_dist_d));
	checkCudaErrors(cudaFree(res_sigma_d));
	checkCudaErrors(cudaFree(res_seq_d));
	checkCudaErrors(cudaFree(dummy_res_dist_out));
	checkCudaErrors(cudaFree(res_seq_out));
	checkCudaErrors(cudaFree(visited_d));
		
	//printf("Outside for loop free memory 2 \n");	
	printf("postprocess time is %f\n",total_time);
	g_active_time += total_time;

}

//(R_d, C_d, d_row, sigma_row, S_row, index, visited_row, local_nodes[tid], left_nodes[tid], right_nodes[tid]);
/*__device__ void check_neighbours(int *&R, int *&C, int *&d_row, int *&dummy_d_row, unsigned long long *&sigma_row, int &node, 
int &left, int &right)
{
		int prv_node = node;
		int a = C[R[node]];
		int b = C[R[node] + 1];
		int dist = 1;
			
		while(a!=left && a!=right)
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
		while(b!=right && b!=left)
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
}*/

__device__ void check_neighbours(int *&R, int *&C, int *&d_row, int *&dummy_d_row, unsigned long long *&sigma_row, int &node,
 int &left, int &right,  bool *&visited)
{
		//int index = 0;
		int a = C[R[node]];
		int b = C[R[node] + 1];
		//S_row[index++] = node;
		visited[node] = true;
		d_row[node] = 0;
		sigma_row[node] = 1;
		int dist = 0;

		while(a!=left && a!=right)
		{
			visited[a] = true;
		//	S_row[index++] = a;
			d_row[a] = ++dist;
			dummy_d_row[a] = d_row[a];
			sigma_row[a] = 1;
			if(visited[C[R[a]]]==false)
			{
				a = C[R[a]];
			}else
			{
				a = C[R[a] + 1];
			}	
		}

		dist = 0;
		while(b!=right && b!=left)
		{
			visited[b] = true;
		//	S_row[index++] = b;
			d_row[b] = ++dist;
			dummy_d_row[b] = d_row[b];
			sigma_row[b] = 1;
			if(visited[C[R[b]]]==false)
			{
				b = C[R[b]];
			}else
			{
				b = C[R[b] + 1];
			}	
		}
}

__global__ void bc_postprocess1(float* bc, unsigned long long *sigma_ld,unsigned long long *sigma_rd, int* d_ld, int* d_rd, int* S_ld,
	int *S_rd, size_t pitch_lsigma, size_t pitch_rsigma, size_t pitch_ld,size_t pitch_rd, size_t pitch_lS, size_t pitch_rS, int *left_dist, 
	int *right_dist, int *local_nodes, int *left_nodes, int *right_nodes, int *R_d, int *C_d,int limit, int n, int *position,int lmemsize, 
	int start_pos_offset, float *delta, size_t pitch_delta, int *res_dist, size_t pitch_res_dist, unsigned long long *res_sigma, 
	size_t pitch_res_sigma, int *dummy_res_dist, int *res_seq, bool *visited,size_t pitch_visited)
{
	int tid = threadIdx.x;
	int index = 0, left, right;
	int *d_lrow, *d_rrow, *S_lrow, *S_rrow;
	unsigned long long *sigma_lrow, *sigma_rrow;
	//if(tid==0)
	//	printf("Outside if tid %d limit %d\n",tid,limit );

	__syncthreads();
	if(blockIdx.x < limit)
	{	
	//	printf("tid %d limit %d\n",tid,limit );
		int *res_dist_row = (int*)((char*)res_dist + blockIdx.x * pitch_res_dist);
		//int *dummy_res_dist_row = (int*)((char*)dummy_res_dist+ blockIdx.x * pitch_dummy_res_dist);
		int *dummy_res_dist_row = (int*)(dummy_res_dist+ blockIdx.x * n);
		//int *res_seq_row = (int*)((char*)res_seq + blockIdx.x * pitch_res_seq);
		int *res_seq_row = (int*)(res_seq + blockIdx.x * n);
		unsigned long long *res_sigma_row = (unsigned long long*)((char*)res_sigma + blockIdx.x * pitch_res_sigma);	
		float *delta_row = (float*)((char*)delta + blockIdx.x * pitch_delta);
		bool *visited_row = (bool*)((char*)visited +  blockIdx.x * pitch_visited);
		
		for(int k=threadIdx.x; k < n; k+=blockDim.x)
		{
			dummy_res_dist_row[k] = INT_MAX;
			//visited_row[k] = false;
			delta_row[k] = 0;	
			res_dist_row[k] = INT_MAX;
			visited_row[k] = false;	
		}
	
		int ind = blockIdx.x;
		left = left_nodes[ind];
		right = right_nodes[ind];

		d_lrow = (int*)((char*)d_ld + (position[left] - start_pos_offset) * pitch_ld);
		sigma_lrow = (unsigned long long*)((char*)sigma_ld + (position[left] - start_pos_offset) * pitch_lsigma);
		S_lrow = (int*)((char*)S_ld + (position[left] - start_pos_offset) * pitch_lS);
	
		if(position[right] < (start_pos_offset + lmemsize))
		{
			d_rrow = (int*)((char*)d_ld + (position[right] - start_pos_offset) * pitch_ld);
			sigma_rrow = (unsigned long long*)((char*)sigma_ld + (position[right] - start_pos_offset) * pitch_lsigma);
			S_rrow = (int*)((char*)S_ld + (position[right] - start_pos_offset) * pitch_lS);	
		}else
		{
			d_rrow = (int*)((char*)d_rd + (position[right] - (start_pos_offset + lmemsize) ) * pitch_rd);
			sigma_rrow = (unsigned long long*)((char*)sigma_rd + (position[right] - (start_pos_offset + lmemsize) ) * pitch_rsigma);
			S_rrow = (int*)((char*)S_rd + (position[right] - (start_pos_offset + lmemsize) ) * pitch_rS);		
		}

		__syncthreads();
	
		for(int k=threadIdx.x; k<n; k+=blockDim.x)
		{
			if((d_lrow[k]+ left_dist[ind]) < (d_rrow[k] + right_dist[ind] ))
			{
				res_dist_row[k] = d_lrow[k]+ left_dist[ind];
				res_sigma_row[k] = sigma_lrow[k];
			}
			else if((d_lrow[k]+ left_dist[ind]) > (d_rrow[k] + right_dist[ind] ))	
			{
				res_dist_row[k] = d_rrow[k] + right_dist[ind];
				res_sigma_row[k] = sigma_rrow[k];
			}
			else
			{
				res_dist_row[k] = d_rrow[k] + right_dist[ind];
				res_sigma_row[k] = sigma_rrow[k] + sigma_lrow[k];
			}
			dummy_res_dist_row[k] = res_dist_row[k];	
			res_seq_row[k] = k;

		//	if(local_nodes[blockIdx.x]==4)
		//		printf("dummy_res_dist_row[%d] = %d , res_seq_row[%d] = %d blockIdx %d\n",k,dummy_res_dist_row[k],k,res_seq_row[k], blockIdx.x);
		}	
	
		__syncthreads();

		if(threadIdx.x==0)
		{
			dummy_res_dist_row[local_nodes[blockIdx.x]] = 0;
			res_dist_row[local_nodes[blockIdx.x]] = 0;	
			res_sigma_row[local_nodes[blockIdx.x]] = 1;
		//	check_neighbours(R_d, C_d, res_dist_row, dummy_res_dist_row, res_sigma_row, local_nodes[blockIdx.x],
		//	left_nodes[blockIdx.x], right_nodes[blockIdx.x], visited_row);
			{
					//int index = 0;
				int a = C_d[R_d[local_nodes[blockIdx.x]]];
				int b = C_d[R_d[local_nodes[blockIdx.x]] + 1];
				//S_row[index++] = node;
				visited_row[local_nodes[blockIdx.x]] = true;
				//res_dist_row[local_nodes[blockIdx.x]] = 0;
				//res_sigma_row[local_nodes[blockIdx.x]] = 1;
			//	printf("first a = %d b = %d loc1 = %d loc2 = %d node %d\n",a,b,R_d[local_nodes[blockIdx.x]],R_d[local_nodes[blockIdx.x]] + 1, local_nodes[blockIdx.x]);
				int dist = 0;
			//	while(a!=left_nodes[blockIdx.x] && a!=right_nodes[blockIdx.x])
				while((R_d[a+1]-R_d[a])==2)
				{
			//		printf("aaa = %d bbb = %d node = %d\n",a,b,local_nodes[blockIdx.x]);
					visited_row[a] = true;
					if(dist > 10)
					{	
						printf("a# = %d left = %d right = %d node = %d\n",a,left_nodes[blockIdx.x],right_nodes[blockIdx.x], local_nodes[blockIdx.x] );
						break;
					}	
				//	S_row[index++] = a;
					res_dist_row[a] = ++dist;
					dummy_res_dist_row[a] = res_dist_row[a];
					res_sigma_row[a] = 1;
					if(visited_row[C_d[R_d[a]]]==false)
					{
						a = C_d[R_d[a]];
					}else
					{
						a = C_d[R_d[a] + 1];
					}	
				}
				dist = 0;
				//while(b!=right_nodes[blockIdx.x] && b!=left_nodes[blockIdx.x])
				while((R_d[b+1]-R_d[b])==2)
				{
			//		printf("aaaaa = %d bbbbb = %d node = %d\n",a,b, local_nodes[blockIdx.x]);
					visited_row[b] = true;
					if(dist > 10)
					{	
						printf("b# = %d left = %d right = %d node = %d\n",b,left_nodes[blockIdx.x],right_nodes[blockIdx.x], local_nodes[blockIdx.x]);
						break;
					}	
				//	S_row[index++] = b;
					res_dist_row[b] = ++dist;
					dummy_res_dist_row[b] = res_dist_row[b];
					res_sigma_row[b] = 1;
					if(visited_row[C_d[R_d[b]]]==false)
					{
						b = C_d[R_d[b]];
					}else
					{
						b = C_d[R_d[b] + 1];
					}	
				}

			}

		//	printf("dummy_res_dist_row[%d] = %d , res_seq_row[%d] = %d tid %d blockid %d\n",local_nodes[blockIdx.x],dummy_res_dist_row[local_nodes[blockIdx.x]],local_nodes[blockIdx.x],res_seq_row[local_nodes[blockIdx.x]],threadIdx.x,blockIdx.x);
		//	printf("node= %d  left %d right %d\n", local_nodes[blockIdx.x], left_nodes[blockIdx.x], right_nodes[blockIdx.x]);
		}	
		__syncthreads();
	}	

}

__global__ void bc_postprocess2(float *bc, int *R, int *C, int g_n, float *delta_d, size_t pitch_delta, int *res_dist, 
unsigned long long *res_sigma, int *res_seq_out, int *dummy_res_dist_out, size_t pitch_res_dist, size_t pitch_res_sigma,
int *offset_Q, size_t pitch_offset_Q, int m, int *F, int *reachval, int limit, float *tmp_bc)
{
	int tid = threadIdx.x;
	//int *dummy_res_dist_out_row = (int *)((char*)dummy_res_dist_out + blockIdx.x * pitch_dummy_res_dist_out);
	int *dummy_res_dist_out_row = (int *)(dummy_res_dist_out + blockIdx.x * g_n);
	//int *res_seq_out_row = (int *)((char*)res_seq_out + blockIdx.x * pitch_res_seq_out );
	if(blockIdx.x < limit)		
	{
		int *d_row = (int*)((char*)res_dist + blockIdx.x*pitch_res_dist);
		unsigned long long *sigma_row = (unsigned long long*)((char*)res_sigma + blockIdx.x*pitch_res_sigma);
		float *delta_row = (float*)((char*)delta_d + blockIdx.x*pitch_delta);
		//int *S_row = (int*)((char*)res_seq_out + blockIdx.x*pitch_res_seq_out);
		int *S_row = (int*)(res_seq_out + blockIdx.x*g_n);
		

		__shared__ int *offset_Q_row;
		__shared__ int offset_ends;
		__shared__ int current_depth;
		__shared__ int start;
	
		if(tid==0)
		{
			offset_Q_row = (int *)((char*)offset_Q + pitch_offset_Q * blockIdx.x);
			offset_Q_row[0] = 0;
			offset_Q_row[ dummy_res_dist_out_row[g_n-1] + 1] = g_n;  
			offset_ends = dummy_res_dist_out_row[g_n-1] + 1;
			current_depth = offset_ends - 1;
			start = S_row[0];
		//	printf("start vertex is %d\n", start);
			//printf("start %d g_n %d seq %d %d %d %d %d %d blockIdx %d current_depth %d\n",start,g_n,S_row[0], S_row[1], S_row[2], S_row[3], S_row[4], S_row[5], blockIdx.x, current_depth);
		}

		__syncthreads();

		for(int i=threadIdx.x; i<g_n; i+=blockDim.x)
		{
			if(i!=0)
			{
				if(dummy_res_dist_out_row[i]!=dummy_res_dist_out_row[i-1])
				{
					offset_Q_row[ dummy_res_dist_out_row[i] ] = i; 
					//if(start==4)
					//	printf("offset_Q_row[%d] val %d\n",dummy_res_dist_out_row[i],i );
				}
			}
		}

		__syncthreads();

		while(current_depth > 0)
		{
			int stack_iter_len = offset_Q_row[current_depth+1]-offset_Q_row[current_depth];
			//if(start==4)
		//	printf("stack_iter_len %d tid %d start %d\n",stack_iter_len,tid, start);
			if(stack_iter_len>512)
			{
				for(int kk=threadIdx.x; kk<2*m; kk+=blockDim.x)
				{
					int w = F[kk];
					if(d_row[w] == current_depth)
					{
						int v = C[kk];
						if(d_row[v] == (d_row[w]+1))
						{
					//	printf("parent %d child %d reachval %d start %d current_depth %d\n",w,v,reachval[w],i, current_depth );
							float change = (sigma_row[w]/(float)sigma_row[v])*(1.0f+delta_row[v]);
							if(reachval[v]!=-1)
							{
								change += ((sigma_row[w]/(float)sigma_row[v])*(1.0f)) * reachval[v] * 2;
							}	
							
							/*if(current_depth==1 && reachval[start]!=-1)
							{
								change += change * reachval[start] * 2;								
							}*/

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
					float sw = (float)sigma_row[w];
					for(int z=R[w]; z<R[w+1]; z++)
					{
						int v = C[z];
						if(d_row[v] == (d_row[w]+1))
						{
							if(reachval[v]!=-1)
							{
								dsw += ((sw/(float)sigma_row[v])*(1.0f)) * reachval[v] * 2;
							//	printf("Inside IFFF parent %d child %d reachval %d start %d current_depth %d sigma[%d] %u dsw %f t_nodes %d\n",w,v,reachval[v],start,current_depth,v,sigma_row[v],dsw,g_n);
							}
							dsw += (sw/(float)sigma_row[v])*(1.0f+delta_row[v]);
						//	printf("parent %d child %d reachval %d start %d current_depth %d sigma[%d] %u dsw %f t_nodes %d\n",w,v,reachval[v],start,current_depth,v,sigma_row[v],dsw,g_n);
						//	if(start==4)
						//		printf("parent %d child %d reachval %d start %d current_depth %d sigma[%d] %u dsw %f t_nodes %d\n",w,v,reachval[v],start,current_depth,v,sigma_row[v],dsw,g_n);
						}
					}
					delta_row[w] += dsw;
					/*if(current_depth==1 && reachval[i]!=-1)
					{
						delta_row[w] += delta_row[w] * reachval[i] * 2;								
					}*/	
				}
			}
			__syncthreads();
			if(tid == 0)
			{
			//	printf("tid %d current_depth %d blockIdx.x %d stack_iter_len %d start %d\n",tid,current_depth,blockIdx.x,stack_iter_len,start );
				current_depth--;
			}
			__syncthreads();
		}

		for(int kk=threadIdx.x; kk<g_n; kk+=blockDim.x)
		{
			atomicAdd(&bc[kk],delta_row[kk]); //Would need to check that kk != i here, but delta_row[kk] is guaranteed to be 0.
			//atomicAdd(&tmp_bc[kk],delta_row[kk]);
		}
		
	}

}

/*
for(int k=threadIdx.x; k < (n*2); k+=blockDim.x)
	{
		if(k < n)
		{
			key_dist_row[k] = d_lrow[k];
			value_nodes_row[k] = k; 
		}else
		{
			key_dist_row[k] = d_rrow[k-n];
			value_nodes_row[k] = k;
		}
	}
	__synchthreads();

	for(int k=threadIdx.x; k < (n*2); k+=blockDim.x)
	{
		atomicMin(&res_dist_row[value_nodes_row[k]],key_dist_row[k]);
	}
*/

