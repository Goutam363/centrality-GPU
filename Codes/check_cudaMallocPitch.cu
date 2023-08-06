#include <stdio.h>
#include <iostream>
#include <cuda.h>

using namespace std;

int main()
{

	int n= 5000;
	int m = 26;
	//int arr[5000][5000];
	int arr_h[n][m];
	int *arr_d;
	int *brr_d;
	size_t pitch_arr, pitch_brr;

	cudaMallocPitch(&arr_d, &pitch_arr, sizeof(int)*n,m);

	dimGrid.x = 13;
	dimGrid.y = 1;
	dimGrid.z = 1;

	dimBlock.x = 1024;
	dimBlock.y = 1;
	dimBlock.z = 1;

	kernel1<<<dimGrid,dimBlock>>>(arr_d,n,m,pitch_arr);
	cudaDeviceSynchronize();

	brr_d = arr_d;
	pitch_brr = pitch_arr;

	kernel2<<<dimGrid, dimBlock>>>(brr_d,n,m,pitch_brr);
	cudaDeviceSynchronize();

	//cudaMemcpy2D(arr_h,sizeof(int)*n, brr_d, pitch_brr, sizeof(int));	
	printf("done\n");

}


__global__ void kernel1(int *a, int n, int m, size_t pitch_arr)
{

	tid = threadId.x;

	if(tid < n )
	{
		for(int itr = 0; itr<14; itr+=13)
		{	
			int *a_row = (int *)((char *)a + (blockId.x + itr) *pitch_arr);
			__syncthreads();
			for(int k = threadId.x; k < n; k+=blockDim.x)
			{
				a_row[k] = 0;
			}	
			__syncthreads();
			for(int k= threadId.x; k < n; k+= blockDim.x)
			{
				a_row[k] += 2;
			}	
		}	
	}	
}

__global__ void kernel2(int *a, int n, int m, size_t pitch_arr)
{

	tid = threadId.x;

	if(tid < n )
	{
		for(int itr = 0; itr<14; itr+=13)
		{	
			int *a_row = (int *)((char *)a + (blockId.x + itr)*pitch_arr);
			__syncthreads();
			for(int k = threadId.x; k < n; k+=blockDim.x)
			{
				a_row[k] += 1;
			}	

			if(threadId.x==0)
				printf("threadId %d blockId %d value %d\n",tid,blockId.x,a_row[threadId.x]);
		}	
	}	
}