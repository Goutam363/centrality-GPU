#include "kernels.cuh"
//#include "ear.h"

#define DIAMETER_SAMPLES 512

void print_reordering(std::vector<int>reordering,std::vector<int>nomap)
{

	printf("Reordring Print\n");
	for(int i=0;i<reordering.size();i++)
	{
		printf("vertex = %d\n",nomap[reordering[i]]);
	}
}

//For portability reasons, we will not use CUDA 6 features here.

void bc_gpu_active(int n, int m, int *&R, int *&C, int *&F, int max_threads_per_block, int number_of_SMs, const std::set<int> &source_vertices,
	std::vector<int>&parents, std::vector<int>&levels, std::vector<int>&level_order, std::vector<int>&level_offset, ear_premble **&ear_info, 
	int max_earlevel, float*&bc_gpu, double&g_active_time, int *&ear_R, int *&ear_C, std::vector<int>&reachval_h)
{
	//Host result data
	
	//float *bc_gpu = new float[g.n];
	float *tmp_bc_gpu = new float[n];

	int *next_source = new int;
	GpuTimer timer;	

	//std::cout<<"Memory allocation has started "<<std::endl;
	printf("max_threads_per_block %d number_of_SMs %d\n",max_threads_per_block,number_of_SMs);
	//Device pointers
	int total_act_vert_size = 0;
	float *bc_d, *delta_d, *tmp_bc_d;
	int *d_d,*d_ld, *R_d, *C_d, *F_d, *Q_d, *Q2_d,*S_d,*S_ld, *endpoints_d, *next_source_d;
	int *d_rd,*S_rd;
	int *sigma_d;
	int *sigma_ld;
	int *sigma_rd;
	size_t pitch_d,pitch_ld, pitch_sigma,pitch_lsigma, pitch_delta, pitch_Q, pitch_Q2,pitch_S,pitch_lS,pitch_endpoints;
	size_t pitch_rd,pitch_rsigma,pitch_rS;

	float G_time = 0;
	double total_time=0;	
	cudaEvent_t sstart,eend;	

	int *jia_d, *diameters_d;
	int left_offset=0,lmemsize=0,rmemsize=0,act_vert_size,reusable_memory,start_pos_offset=0,position_index = 0;
	//Grid parameters
	dim3 dimBlock, dimGrid;
	dimGrid.x = number_of_SMs;
	dimGrid.y = 1;
	dimGrid.z = 1;

	dimBlock.x = max_threads_per_block;
	dimBlock.y = 1;
	dimBlock.z = 1;

	fflush(stdout);
	int *dummy = (int *)malloc(sizeof(int) * 5);
	fflush(stdout);
	std::vector<int>position(n,0);
	next_source[0] = number_of_SMs; 

//	arrange_by_threshold(level_offset,level_order);
	//Allocate and transfer data to the GPU
	thrust::device_vector<int> reachval_d = reachval_h;

	checkCudaErrors(cudaMalloc((void**)&bc_d,sizeof(float)*n));
	//checkCudaErrors(cudaMalloc((void**)&tmp_bc_d,sizeof(float)*n));
	checkCudaErrors(cudaMalloc((void**)&R_d,sizeof(int)*(n+1)));
	checkCudaErrors(cudaMalloc((void**)&C_d,sizeof(int)*(2*m)));
	checkCudaErrors(cudaMalloc((void**)&F_d,sizeof(int)*(2*m)));

	checkCudaErrors(cudaMemcpy(R_d,R,sizeof(int)*(n+1),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(C_d,C,sizeof(int)*(2*m),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(F_d,F,sizeof(int)*(2*m),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(bc_d,bc_gpu,sizeof(float)*n,cudaMemcpyHostToDevice));
	
	checkCudaErrors(cudaMallocPitch((void**)&d_d,&pitch_d,sizeof(int)*n,dimGrid.x));
	checkCudaErrors(cudaMallocPitch((void**)&sigma_d,&pitch_sigma,sizeof(int)*n,dimGrid.x));
	checkCudaErrors(cudaMallocPitch((void**)&S_d,&pitch_S,sizeof(int)*n,dimGrid.x));

	checkCudaErrors(cudaMallocPitch((void**)&delta_d,&pitch_delta,sizeof(float)*n,dimGrid.x));
	checkCudaErrors(cudaMallocPitch((void**)&Q_d,&pitch_Q,sizeof(int)*n,dimGrid.x)); //Making Queues/Stack of size O(n) since we won't duplicate
	checkCudaErrors(cudaMallocPitch((void**)&Q2_d,&pitch_Q2,sizeof(int)*n,dimGrid.x));
	checkCudaErrors(cudaMallocPitch((void**)&endpoints_d,&pitch_endpoints,sizeof(int)*(n+1),dimGrid.x));

	checkCudaErrors(cudaMalloc((void**)&next_source_d,sizeof(int)));
	std::vector<int>reordering_map;
	checkCudaErrors(cudaMalloc((void**)&jia_d,sizeof(int)));
	checkCudaErrors(cudaMalloc((void**)&diameters_d,sizeof(int)*DIAMETER_SAMPLES));
	thrust::device_vector<int> active_vertices;

	double postproc_time = 0;
	int prv_size=0;

	for(int lev=1;lev<level_offset.size();lev++)
	{
		rmemsize = 0;
	//	checkCudaErrors(cudaMemset(tmp_bc_d,0,sizeof(float)*n));
		checkCudaErrors(cudaMemcpy(next_source_d,next_source,sizeof(int),cudaMemcpyHostToDevice));	
		
		act_vert_size = level_offset[lev]-level_offset[lev-1];
		std::vector<int>reordering(level_order.begin()+level_offset[lev-1],level_order.begin()+level_offset[lev-1]+act_vert_size);
		
		for(int i=0;i<reordering.size();i++)
		{
			position[ reordering[i] ] = position_index + i;//start_pos_offset + i;
		}	

		rmemsize = act_vert_size;// - lmemsize;
		position_index += rmemsize;
		active_vertices = reordering;
		
		total_act_vert_size += act_vert_size;
		//rmemsize = 50;		//To check speedup
		//exit(0);	
		checkCudaErrors(cudaMallocPitch((void**)&d_rd,&pitch_rd,sizeof(int)*n,rmemsize));
		checkCudaErrors(cudaMallocPitch((void**)&sigma_rd,&pitch_rsigma,sizeof(int)*n,rmemsize));
		checkCudaErrors(cudaMallocPitch((void**)&S_rd,&pitch_rS,sizeof(int)*n,rmemsize));

		checkCudaErrors(cudaMemset(jia_d,0,sizeof(int)));
		checkCudaErrors(cudaMemset(diameters_d,0,sizeof(int)*DIAMETER_SAMPLES));
		
		//start_clock(sstart,eend);
		timer.Start();			
		bc_kgpu_active<<<dimGrid,dimBlock>>>(bc_d,R_d,C_d,F_d,n,m,d_d,d_ld,d_rd,sigma_d,sigma_ld,sigma_rd,delta_d,Q_d,Q2_d,S_d,
		S_ld,S_rd,endpoints_d,next_source_d,pitch_d,pitch_ld,pitch_rd,pitch_sigma,pitch_lsigma,pitch_rsigma,pitch_delta,pitch_Q,
		pitch_Q2,pitch_S,pitch_lS,pitch_rS,pitch_endpoints,0,act_vert_size,jia_d,diameters_d, 
		thrust::raw_pointer_cast(active_vertices.data()),rmemsize, thrust::raw_pointer_cast(reachval_d.data()));

		checkCudaErrors(cudaPeekAtLastError());
		//checkCudaErrors(cudaMemcpy(tmp_bc_gpu,tmp_bc_d,sizeof(float)*n,cudaMemcpyDeviceToHost));
		//checkCudaErrors(cudaDeviceSynchronize());
		//G_time = end_clock(sstart,eend);
		timer.Stop();
		total_time += timer.Elapsed();	
			
		//g_active_time += timer.Elapsed();
		//printf("outside kernel.cu test1 %d\n",lev);
		//total_time += G_time;
//		printf("G_time is %f\n",total_time);
		
		active_vertices.clear();
		active_vertices.shrink_to_fit();
	
		//exit(0);				
		if(lev!=1)
		{	

			//printf("Going for postprocess \n");
			postprocess(n, m, R_d, C_d, F_d, bc_d, d_ld, d_rd, sigma_ld, sigma_rd, S_ld, S_rd, lmemsize, ear_info, number_of_SMs, 
			max_threads_per_block, pitch_ld, pitch_rd, pitch_lsigma, pitch_rsigma, pitch_lS, pitch_rS, start_pos_offset, level_order,
			rmemsize, position, postproc_time, ear_R, ear_C, reachval_d, g_active_time);

			//printf("After postprocess \n");
		}
		
		start_pos_offset += lmemsize;

		if(lev!=1)
		{
			checkCudaErrors(cudaFree(d_ld));
			checkCudaErrors(cudaFree(sigma_ld));
			checkCudaErrors(cudaFree(S_ld));
		}
			
		sigma_ld = sigma_rd;
		d_ld = d_rd;
		S_ld = S_rd;

		lmemsize = rmemsize;
		
		pitch_lsigma = pitch_rsigma;
		pitch_ld = pitch_rd;
		pitch_lS = pitch_rS;
	//Transfer result to CPU
	}

	checkCudaErrors(cudaMemcpy(bc_gpu,bc_d,sizeof(float)*n,cudaMemcpyDeviceToHost));

	/*printf("After LOOP PArtial betweenness centrality kernels_active \n");
	for(int h=0; h<n; h++)
	{
		printf("%d = %f\n",h, bc_gpu[h]);
	}*/

//	printf("CUDA FREE\n");
	//Free memory
	checkCudaErrors(cudaFree(bc_d));
	//checkCudaErrors(cudaFree(tmp_bc_d));
	checkCudaErrors(cudaFree(R_d));
	checkCudaErrors(cudaFree(C_d));
	checkCudaErrors(cudaFree(F_d));
	checkCudaErrors(cudaFree(d_d));
	checkCudaErrors(cudaFree(sigma_d));
	checkCudaErrors(cudaFree(delta_d));
	checkCudaErrors(cudaFree(Q_d));
	checkCudaErrors(cudaFree(Q2_d));
	checkCudaErrors(cudaFree(S_d));
	checkCudaErrors(cudaFree(endpoints_d));
	checkCudaErrors(cudaFree(next_source_d));

	checkCudaErrors(cudaFree(d_rd));
	checkCudaErrors(cudaFree(sigma_rd));
	checkCudaErrors(cudaFree(S_rd));

	checkCudaErrors(cudaFree(jia_d));
	checkCudaErrors(cudaFree(diameters_d));

	reachval_d.clear();
	reachval_d.shrink_to_fit();	

	//delete[] bc_gpu;
	delete next_source;
	printf("total time for active vertices is %f\n",total_time);
	g_active_time += total_time;
	//g_active_time += postproc_time;
	//printf("post proc time is %f\n",postproc_time);
	printf("total_act_vert_size is %d\n",total_act_vert_size);
	//return bc_agpu_v;
}

//Note: N must be a power of two
//Simple/Naive bitonic sort. We're only sorting ~512 elements one time, so performance isn't important
__device__ void bitonic_sort_active(int *values, int N)
{
	unsigned int idx = threadIdx.x;

	for (int k = 2; k <= N; k <<= 1)
	{
		for (int j = k >> 1; j > 0; j = j >> 1)
		{
			while(idx < N) 
			{
				int ixj = idx^j;
				if (ixj > idx) 
				{
					if ((idx&k) == 0 && values[idx] > values[ixj]) 
					{
						//exchange(idx, ixj);
						int tmp = values[idx];
						values[idx] = values[ixj];
						values[ixj] = tmp;
					}
					if ((idx&k) != 0 && values[idx] < values[ixj]) 
					{
						//exchange(idx, ixj);
						int tmp = values[idx];
						values[idx] = values[ixj];
						values[ixj] = tmp;
					}
				}
				idx += blockDim.x;
			}
			__syncthreads();
			idx = threadIdx.x;
		}
	}
}

__global__ void bc_kgpu_active(float *bc, const int *R, const int *C, const int *F, const int n, const int m,int*d, int *ld,
int *rd,int *sigma,int *lsigma,int *rsigma, float *delta, int *Q, int *Q2,int*S,
int *lS,int *rS,int *endpoints, int *next_source,size_t pitch_d, size_t pitch_ld,size_t pitch_rd,size_t pitch_sigma,
size_t pitch_lsigma,size_t pitch_rsigma,size_t pitch_delta, size_t pitch_Q, size_t pitch_Q2,size_t pitch_S,size_t pitch_lS,
size_t pitch_rS,size_t pitch_endpoints, int start, int end, int *jia, int *diameters,
int *active_vertices,int rmemsize, int *reachval)
{
	__shared__ int ind ;
	__shared__ int i;
	int j = threadIdx.x;
	int *d_row = (int*)((char*)d + blockIdx.x*pitch_d);
	int *sigma_row = (int*)((char*)sigma + blockIdx.x*pitch_sigma);
	float *delta_row = (float*)((char*)delta + blockIdx.x*pitch_delta);

	int *rd_row,*rS_row;
	int *rsigma_row;

	__shared__ int *Q_row;
	__shared__ int *Q2_row;
	__shared__ int *S_row;
	__shared__ int *endpoints_row;

	if(j == 0)
	{
		{
			ind = blockIdx.x + start;
			//i = ind;
			if(ind < rmemsize)
			{	
				i = active_vertices[ind];
				Q_row = (int*)((char*)Q + blockIdx.x*pitch_Q);
				Q2_row = (int*)((char*)Q2 + blockIdx.x*pitch_Q2);
				S_row = (int*)((char*)S + blockIdx.x*pitch_S);
				endpoints_row = (int*)((char*)endpoints + blockIdx.x*pitch_endpoints);
				*jia = 0;
			}	
		}
	}
	__syncthreads();

	if(ind < rmemsize)
	{
		rS_row = (int*)((char*)rS + ind * pitch_rS);
		rd_row = (int*)((char*)rd + ind * pitch_rd);
		rsigma_row = (int *)((char*)rsigma + ind * pitch_rsigma);	
	}	
	//__syncthreads();
	if((ind==0) && (j < DIAMETER_SAMPLES))
	{
		diameters[j] = INT_MAX;
	}
	__syncthreads();

	while(ind < end)
	{
		for(int k=threadIdx.x; k<n; k+=blockDim.x)
		{
			if(k == i) //If k is the source node...
			{
				d_row[k] = 0;
				sigma_row[k] = 1;
			}
			else
			{
				d_row[k] = INT_MAX;
				sigma_row[k] = 0;
			}	
			delta_row[k] = 0;
		}
		__syncthreads();

		//Shortest Path Calculation
		__shared__ int Q_len;
		__shared__ int Q2_len;
	        __shared__ int S_len;
	        __shared__ int current_depth; 
		__shared__ int endpoints_len;
		__shared__ bool sp_calc_done;

		if(j == 0)
		{
			Q_row[0] = i;
			Q_len = 1;
			Q2_len = 0;
			S_row[0] = i;
			S_len = 1;
			endpoints_row[0] = 0;
			endpoints_row[1] = 1;
			endpoints_len = 2;
			current_depth = 0;
			sp_calc_done = false;
		}
		__syncthreads();

		//Do first iteration separately since we already know the edges to traverse
		for(int r=threadIdx.x+R[i]; r<R[i+1]; r+=blockDim.x)
		{
			int w = C[r];
			//No multiple/self edges - each value of w is unique, so no need for atomics
			if(d_row[w] == INT_MAX)
			{
				d_row[w] = 1; 
				int t = atomicAdd(&Q2_len,1);
				Q2_row[t] = w;
			}
			if(d_row[w] == (d_row[i]+1))
			{
				atomicAdd(&sigma_row[w],1); 
			}
		}
		__syncthreads();

		if(Q2_len == 0)
		{
			sp_calc_done = true;
		}
		else
		{
			for(int kk=threadIdx.x; kk<Q2_len; kk+=blockDim.x)
			{
				Q_row[kk] = Q2_row[kk];
				S_row[kk+S_len] = Q2_row[kk];
			}
			__syncthreads();
			if(j == 0)
			{
				endpoints_row[endpoints_len] = endpoints_row[endpoints_len-1] + Q2_len;
				endpoints_len++;
				Q_len = Q2_len;
				S_len += Q2_len;
				Q2_len = 0;
				current_depth++;
			}
		}
		__syncthreads();

		while(!sp_calc_done)
		{
			if((*jia) && (Q_len > 512))
			{
				for(int k=threadIdx.x; k<2*m; k+=blockDim.x)
				{
					int v = F[k];
					if(d_row[v] == current_depth) 
					{
						int w = C[k];
						if(atomicCAS(&d_row[w],INT_MAX,d_row[v]+1) == INT_MAX)
						{
							int t = atomicAdd(&Q2_len,1);
							Q2_row[t] = w;
						}
						if(d_row[w] == (d_row[v]+1))
						{
							atomicAdd(&sigma_row[w],sigma_row[v]);
						}
					}	
				}
			}
			else
			{
				__shared__ int next_index;
				if(j == 0)
				{
					next_index = blockDim.x;
				}
				__syncthreads();
				int k = threadIdx.x; //Initial vertices
				while(k < Q_len)
				{
					int v = Q_row[k];
					for(int r=R[v]; r<R[v+1]; r++)
					{
						int w = C[r];
						//Use atomicCAS to prevent duplicates
						if(atomicCAS(&d_row[w],INT_MAX,d_row[v]+1) == INT_MAX)
						{
							int t = atomicAdd(&Q2_len,1);
							Q2_row[t] = w;
						}
						if(d_row[w] == (d_row[v]+1))
						{
							atomicAdd(&sigma_row[w],sigma_row[v]);
						}
					}
					k = atomicAdd(&next_index,1);
				}
			}
			__syncthreads();

			if(Q2_len == 0) //If there is no additional work found, we're done
			{
				break;
			}
			else //If there is additional work, transfer elements from Q2 to Q, reset lengths, and add vertices to the stack
			{
				for(int kk=threadIdx.x; kk<Q2_len; kk+=blockDim.x)
				{
					Q_row[kk] = Q2_row[kk];
					S_row[kk+S_len] = Q2_row[kk];
				}
				__syncthreads();
				if(j == 0)
				{
					endpoints_row[endpoints_len] = endpoints_row[endpoints_len-1] + Q2_len;
					endpoints_len++;
					Q_len = Q2_len;
					S_len += Q2_len;
					Q2_len = 0;
					current_depth++;
				}
				__syncthreads();
			}
		}

		//The elements at the end of the stack will have the largest distance from the source
		//Using the successor method, we can start from one depth earlier
		if(j == 0)
		{
			current_depth = d_row[S_row[S_len-1]] - 1;
			if(ind<DIAMETER_SAMPLES)
			{
				diameters[ind] = current_depth+1;
			}
		}
		__syncthreads();

		//Dependency Accumulation (Madduri/Ediger successor method)
		while(current_depth > 0)
		{
			int stack_iter_len = endpoints_row[current_depth+1]-endpoints_row[current_depth];
			if((*jia) && (stack_iter_len>512))
			{
				for(int kk=threadIdx.x; kk<2*m; kk+=blockDim.x)
				{
					int w = F[kk];
					if(d_row[w] == current_depth)
					{
						int v = C[kk];
						if(d_row[v] == (d_row[w]+1))
						{
							float change = (sigma_row[w]/(float)sigma_row[v])*(1.0f * (reachval[v]+1) + delta_row[v]);	
							atomicAdd(&delta_row[w],change);
						}		
					}
				}
			}
			else 
			{
				for(int kk=threadIdx.x+endpoints_row[current_depth]; kk<endpoints_row[current_depth+1]; kk+=blockDim.x)
				{
					int w = S_row[kk];
					float dsw = 0;
					float sw = (float)sigma_row[w];
					for(int z=R[w]; z<R[w+1]; z++)
					{
						int v = C[z];
						if(d_row[v] == (d_row[w]+1))
						{
							dsw += (sw/(float)sigma_row[v]) *  (1.0f * (reachval[v]+1) + delta_row[v]);
							//if(i==1)
							//	printf("start node = %d ppparent %d child %d reachval %d start %d current_depth %d sigma[%d] %u dsw %f t_nodes %d delta_row[%d] = %f sw = %f res = %f partres = %f %d %f %f\n",i,w,v,reachval[v],start,current_depth,v,sigma_row[v],dsw,n,v,delta_row[v],sw,((sw/(float)sigma_row[v])*(1.0f * (reachval[v]+1) + delta_row[v])), (1.0f * (reachval[v]+1) + delta_row[v]),reachval[v],delta_row[v],(1.0f * (reachval[v]+1)));
						}	
					}
					delta_row[w] += dsw;	//When art vertex is start node of spanning tree
				}
			}
			__syncthreads();
			if(j == 0)
			{
				current_depth--;
			}
			__syncthreads();
		}

		for(int kk=threadIdx.x; kk<n; kk+=blockDim.x)
		{
			atomicAdd(&bc[kk],(delta_row[kk]*(reachval[i]+1))); //Would need to check that kk != i here, but delta_row[kk] is guaranteed to be 0.
			//atomicAdd(&tmp_bc[kk],(delta_row[kk]*(reachval[i]+1)));
		}
		
		
		if(ind < rmemsize)
		{
			for(int kk=threadIdx.x; kk<n; kk+=blockDim.x)
			{
					rd_row[kk] 		= 	d_row[kk];
	 				rsigma_row[kk]	=   sigma_row[kk];
 					rS_row[kk]		= 	S_row[kk];
 			}	
		}	

		if(j == 0)
		{
			ind = atomicAdd(next_source,1);

			if(ind < rmemsize)
			{
				{	
						i = active_vertices[ind];
				}
			}		
		}

		__syncthreads();


		if(ind < rmemsize)
		{
			rS_row = (int*)((char*)rS + ind * pitch_rS);
			rd_row = (int*)((char*)rd + ind * pitch_rd);
			rsigma_row = (int *)((char*)rsigma + ind * pitch_rsigma);	
		}	
		
		if(ind == 2*DIAMETER_SAMPLES) //Might want to play around with this number. Safe to assume that they are done by now? Probably...
		{
			__shared__ int diameter_keys[DIAMETER_SAMPLES];
			for(int kk = threadIdx.x; kk<DIAMETER_SAMPLES; kk+=blockDim.x)
			{
				diameter_keys[kk] = diameters[kk];
			}
			__syncthreads();
			bitonic_sort_active(diameter_keys,DIAMETER_SAMPLES);
			__syncthreads();
			if(j == 0)
			{
				int log2n = 0;
				int tempn = n;
				while(tempn >>= 1)
				{
					++log2n;
				}
				if(diameter_keys[DIAMETER_SAMPLES/2] < 4*log2n) //Use the median
				{
					*jia = 1;
				}
			}
		}
		__syncthreads();
	}
}

