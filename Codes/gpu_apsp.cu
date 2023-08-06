/**
 * Modified from orignal version of Pawan Harish .
 *
 * GPU implementation of APSP using multiple SSSP. 
 * Each SSSP is parallel, multipe calls are sequential.
 **/

#include "gpu_apsp.h"
#include "utils.h"

#define MAX_THREADS_PER_BLOCK 512
#define MAX_COST 10000000

#include "DijkastraKernel2.cu"
#include "DijkastraKernel.cu"
#include "BetKernel.cu"


int iteration=0;
int stream_id=0;
int iter=0;
int flag_cpy=0;

void GPU_APSP::run_dummy_memcpy() {
	int n = 1000;
	int *h_array = (int*) malloc(sizeof(int)*n);
	for ( int i = 0 ; i < n; i++ )
		h_array[i] = i;

	int d_n;
	int *d_array;
	if( cudaMalloc( (void**) &d_n, sizeof(int)*n)!=cudaSuccess)
		cout<<"error d_n"<<endl;
	( cudaMemcpy( d_array, h_array, sizeof(int)*n, cudaMemcpyHostToDevice) );
}

void GPU_APSP::run_apsp(int start, int num,int totalvertices,int bccnum) {
	// setup execution parameters
	BiCC bicc;
	dim3  grid( num_of_blocks, 1, 1);
	dim3  threads( num_of_threads_per_block, 1, 1);
	
	APSP_Graph ap;
	GpuTimer timer;
	int k=0;
		
	h_betweenness=(double *)malloc(sizeof(double)*totalvertices);
	memset(h_betweenness,0,sizeof(double)*totalvertices);
	
	/*for(int b=0;b<totalvertices;b++)	
	{
		h_betweenness[b]=ap.betweenness[b];
		cout<<"bet "<<b<<" ="<<h_betweenness[b]<<endl;
	}*/

	if(cudaMalloc((void **)&d_betweenness,sizeof(double)*totalvertices)!=cudaSuccess)
			cout<<"error d_betweenness";
	cudaMemcpy(d_betweenness,h_betweenness,sizeof(double)*totalvertices,cudaMemcpyHostToDevice);
	timer.Start();
	{
		
	
			BetCKernel<<< grid, threads>>>( d_graph_nodes, d_graph_edges, d_graph_weights, no_of_nodes, edge_list_size,start,num_of_blocks,num_of_threads_per_block,d_old_vertex,d_reachval,d_betweenness,d_arr1,d_sigma,d_dist,d_upq,d_vpq,d_pq_index,d_bstack,d_bstackind,d_pindex,d_delta,d_parent_bool);
		
		cudaDeviceSynchronize();
		cudaMemcpy(h_betweenness, d_betweenness, sizeof(double)*totalvertices,cudaMemcpyDeviceToHost);
	//	cout<<"gpu_output"<<endl;
		for(int i=0;i<no_of_nodes;i++)
		{
		//	if(bicc.bccvertno[bccnum][i].second==-1)
				ap.betweenness[h_old_vertex[i]]+=h_betweenness[h_old_vertex[i]];
		//	else
		//		ap.betweenness[h_old_vertex[i]]=h_betweenness[h_old_vertex[i]];			

			//cout<<h_old_vertex[i]<<" = "<<h_betweenness[h_old_vertex[i]]<<endl;
		}
		
		
	}
	timer.Stop();
	printf( "\nGPU Processing time: %f (ms)\n", timer.Elapsed());
}

void GPU_APSP::set_graph (Graph &graph,int bccnum) {
	
	BiCC bicc;
	no_of_nodes = graph.n_vertices;
	edge_list_size = graph.n_edges * 2;

	h_graph_nodes = (int*) malloc(sizeof(int)*no_of_nodes);
	h_graph_edges = (int*) malloc(sizeof(int)*edge_list_size);
	h_upq = (int*) malloc(sizeof(int)*edge_list_size);
	h_vpq = (int*) malloc(sizeof(int)*edge_list_size);
	h_graph_weights = (short int*) malloc(sizeof(short int)*edge_list_size);
	h_reachval = (int*) malloc(sizeof(int)*no_of_nodes);
	h_old_vertex = (int*) malloc(sizeof(int)*no_of_nodes);
	
	int edge_no = 0;
	for ( int i = 0 ; i < no_of_nodes; i++ ) {
		h_graph_nodes[i] = edge_no;
		for ( int j = 0; j < graph.adj_list[i].size(); j++ ) {
			h_graph_edges[edge_no] = graph.adj_list[i][j].first;
			h_graph_weights[edge_no] = graph.adj_list[i][j].second;
			h_upq[edge_no]=i;
			h_vpq[edge_no]=graph.adj_list[i][j].first;
			edge_no++;
		}
	}

	int reach=0;
	//cout<<"gp reach values are"<<endl;
	for(int i=0;i<no_of_nodes;i++)
	{
		h_old_vertex[i]=bicc.bccvertno[bccnum][i].first;
		if(bicc.bccvertno[bccnum][i].second==-1)
			reach=0;
		else
			reach=bicc.bccvertno[bccnum][i].second;
		h_reachval[i]=reach;
		//cout<<"new_vertex= "<<i<<" old_vertex= "<<h_old_vertex[i]<<" reachval= "<<reach<<endl;
			
	}
	num_of_blocks = 1;
	num_of_threads_per_block = no_of_nodes;

	//Make execution Parameters according to the number of nodes
	//Distribute threads across multiple Blocks if necessary
	if(no_of_nodes>MAX_THREADS_PER_BLOCK) {
		num_of_blocks = (int)ceil(no_of_nodes/(double)MAX_THREADS_PER_BLOCK); 
		num_of_threads_per_block = MAX_THREADS_PER_BLOCK; 
	}

	h_graph_mask = (bool*) malloc(sizeof(bool)*no_of_nodes);
	h_graph_updating_cost = (int*) malloc(sizeof(int)*no_of_nodes);

	// initalize the memory
	for( unsigned int i = 0; i < no_of_nodes; i++) {
		h_graph_updating_cost[i] = MAX_COST;
		h_graph_mask[i]=false;
	}

	//Copy h_old_vertrex to d_old_vertex	
	if(cudaMalloc((void**)&d_old_vertex,sizeof(int)*no_of_nodes)!=cudaSuccess)
		cout<<"Error d_old_vertex"<<endl;

	cudaMemcpy(d_old_vertex,h_old_vertex,sizeof(int)*no_of_nodes,cudaMemcpyHostToDevice);

	//Copy h_reachval to d_reachval	
	if(cudaMalloc((void**)&d_reachval,sizeof(int)*no_of_nodes)!=cudaSuccess)
		cout<<"Error d_reachval"<<endl;

	cudaMemcpy(d_reachval,h_reachval,sizeof(int)*no_of_nodes,cudaMemcpyHostToDevice);
	
	
	//Copy the int list to device memory
	if( cudaMalloc( (void**) &d_graph_nodes, sizeof(int)*no_of_nodes)!=cudaSuccess)
		cout<<"Error d_graph_nodes"<<endl;

	( cudaMemcpy( d_graph_nodes, h_graph_nodes, sizeof(int)*no_of_nodes, cudaMemcpyHostToDevice) );

	//Copy the Edge List to device Memory
	if( cudaMalloc( (void**) &d_graph_edges, sizeof(int)*edge_list_size)!=cudaSuccess)
		cout<<"Error d_graph_edges"<<endl;

	( cudaMemcpy( d_graph_edges, h_graph_edges, sizeof(int)*edge_list_size, cudaMemcpyHostToDevice) );

	if( cudaMalloc( (void**) &d_graph_weights, sizeof(short int)*edge_list_size)!=cudaSuccess)
			cout<<"Error d_graph_weights"<<endl;
	( cudaMemcpy( d_graph_weights, h_graph_weights, sizeof(short int)*edge_list_size, cudaMemcpyHostToDevice) );

	int size_tb=num_of_threads_per_block*num_of_blocks;
	
	int limit_nodes=size_tb*no_of_nodes;
	int limit_edges=size_tb*edge_list_size;

	h_arr1=(int*)malloc(sizeof(int)*size_tb*no_of_nodes);	
	for(int e=0;e<limit_nodes;e++)
	{
		h_arr1[e]=NOT_REACHABLE;
	}

	if(cudaMalloc((void**) &d_arr1,sizeof(int)*size_tb*no_of_nodes)!=cudaSuccess)
		cout<<"d_arr1 error"<<endl;
	cudaMemcpy(d_arr1,h_arr1,sizeof(int)*size_tb*no_of_nodes,cudaMemcpyHostToDevice);

	//**arr1=(int **)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);
	//arr1[tid]=(int *)malloc(sizeof(int)*no_of_nodes);

	h_sigma=(int *)malloc(sizeof(int)*size_tb*no_of_nodes);
	memset(h_sigma,0,sizeof(int)*size_tb*no_of_nodes);	
	if(cudaMalloc((void**) &d_sigma,sizeof(int)*size_tb*no_of_nodes)!=cudaSuccess)
		cout<<"d_sigma error"<<endl;
	cudaMemcpy(d_sigma,h_sigma,sizeof(int)*size_tb*no_of_nodes,cudaMemcpyHostToDevice);
	
	//int **sigma=(int **)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);
	//sigma[tid]=(int *)malloc(sizeof(int)*no_of_nodes);
	
	h_bstack=(int *)malloc(sizeof(int)*size_tb*no_of_nodes);
	memset(h_bstack,0,sizeof(int)*size_tb*no_of_nodes);
	if(cudaMalloc((void**) &d_bstack,sizeof(int)*size_tb*no_of_nodes)!=cudaSuccess)
		cout<<"d_bstack error"<<endl;
	cudaMemcpy(d_bstack,h_bstack,sizeof(int)*size_tb*no_of_nodes,cudaMemcpyHostToDevice);

	//int **bstack=(int **)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);
	//bstack[tid]=(int *)malloc(sizeof(int)*no_of_nodes);

	h_bstackind=(int *)malloc(sizeof(int)*size_tb);
	memset(h_bstackind,0,sizeof(int)*size_tb);	
	if(cudaMalloc((void**) &d_bstackind,sizeof(int)*size_tb)!=cudaSuccess)
		cout<<"d_bstackind error"<<endl;
	cudaMemcpy(d_bstackind,h_bstackind,sizeof(int)*size_tb,cudaMemcpyHostToDevice);	
	//int *bstackind=(int *)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);

	h_pindex=(int *)malloc(sizeof(int)*size_tb*no_of_nodes);
	memset(h_pindex,0,sizeof(int)*size_tb*no_of_nodes);
	if(cudaMalloc((void**) &d_pindex,sizeof(int)*size_tb*no_of_nodes)!=cudaSuccess)
		cout<<"d_pindex error"<<endl;
	cudaMemcpy(d_pindex,h_pindex,sizeof(int)*size_tb*no_of_nodes,cudaMemcpyHostToDevice);


//	int **pindex=(int **)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);
//	pindex[tid]=(int *)malloc(sizeof(int)*no_of_nodes);
	
	h_delta=(double *)malloc(sizeof(double)*size_tb*no_of_nodes);
	memset(h_delta,0,sizeof(double)*size_tb*no_of_nodes);
	if(cudaMalloc((void**) &d_delta,sizeof(double)*size_tb*no_of_nodes)!=cudaSuccess)
		cout<<"d_pindex error"<<endl;
	cudaMemcpy(d_delta,h_delta,sizeof(double)*size_tb*no_of_nodes,cudaMemcpyHostToDevice);


	h_dist=(int *)malloc(sizeof(int)*size_tb*edge_list_size);
	for(int e=0;e<limit_edges;e++)
	{
		h_dist[e]=NOT_REACHABLE;
	}	
	if(cudaMalloc((void**) &d_dist,sizeof(int)*size_tb*edge_list_size)!=cudaSuccess)
		cout<<"d_dist error"<<endl;
	cudaMemcpy(d_dist,h_dist,sizeof(int)*size_tb*edge_list_size,cudaMemcpyHostToDevice);
	
	//int **dist=(int **)malloc(num_of_threads_per_block*num_of_blocks*sizeof(int));	
	//dist[tid]=(int *)malloc(edge_list_size*sizeof(int));

	h_upq=(int *)malloc(sizeof(int)*size_tb*edge_list_size);
	memset(h_upq,0,sizeof(int)*size_tb*edge_list_size);
	if(cudaMalloc((void**) &d_upq,sizeof(int)*size_tb*edge_list_size)!=cudaSuccess)
		cout<<"d_upq error"<<endl;
	cudaMemcpy(d_upq,h_upq,sizeof(int)*size_tb*edge_list_size,cudaMemcpyHostToDevice);

	//	int **upq=(int **)malloc(num_of_threads_per_block*num_of_blocks*sizeof(int));	
	//	upq[tid]=(int *)malloc(edge_list_size*sizeof(int));

	h_vpq=(int *)malloc(sizeof(int)*size_tb*edge_list_size);
	memset(h_vpq,0,sizeof(int)*size_tb*edge_list_size);
	if(cudaMalloc((void**) &d_vpq,sizeof(int)*size_tb*edge_list_size)!=cudaSuccess)
		cout<<"d_vpq error"<<endl;
	cudaMemcpy(d_vpq,h_vpq,sizeof(int)*size_tb*edge_list_size,cudaMemcpyHostToDevice);

	//	int **vpq=(int **)malloc(num_of_threads_per_block*num_of_blocks*sizeof(int));	
	//	vpq[tid]=(int *)malloc(edge_list_size*sizeof(int));
	
	h_pq_index=(int *)malloc(sizeof(int)*size_tb);
	for(int e=0;e<size_tb;e++)
	{
		h_pq_index[e]=1;
	}
	if(cudaMalloc((void**) &d_pq_index,sizeof(int)*size_tb)!=cudaSuccess)
		cout<<"d_pq_index error"<<endl;
	cudaMemcpy(d_pq_index,h_pq_index,sizeof(int)*size_tb,cudaMemcpyHostToDevice);

	//	int *pq_index=(int *)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);
	
	/*int *h_t_parent=(int *)malloc(sizeof(int)*size_tb*no_of_nodes*no_of_nodes);
	memset(h_t_parent,0,sizeof(int)*size_tb*no_of_nodes*no_of_nodes);
	if(cudaMalloc((void**) &d_t_parent,sizeof(int)*size_tb*no_of_nodes*no_of_nodes)!=cudaSuccess)
		cout<<"d_t_parent error"<<endl;
	cudaMemcpy(d_t_parent,h_t_parent,sizeof(int)*size_tb*no_of_nodes*no_of_nodes,cudaMemcpyHostToDevice);		*/

	//int ***t_parent=(int ***)malloc(sizeof(int)*num_of_threads_per_block*num_of_blocks);
	//	t_parent[tid]=(int **)malloc(sizeof(int)*no_of_nodes*no_of_nodes);
	//	t_parent[tid][tid]=(int *)malloc(sizeof(int)*no_of_nodes);
	bool *h_parent_bool = (bool*) malloc(sizeof(bool)*edge_list_size*size_tb);
	memset(h_parent_bool,0,sizeof(bool)*size_tb*edge_list_size);
	if(cudaMalloc((void**) &d_parent_bool,sizeof(int)*size_tb*edge_list_size)!=cudaSuccess)
		cout<<"d_parent_bool error"<<endl;
	cudaMemcpy(d_parent_bool,h_parent_bool,sizeof(int)*size_tb*edge_list_size,cudaMemcpyHostToDevice);	

	
}

void GPU_APSP::free_graph(){
	// cleanup memory
	/*free( h_graph_nodes);
	free( h_graph_edges);
	free( h_graph_mask);
	free( h_graph_weights);
	free( h_graph_updating_cost);
	free( h_cost);
	free (dist);
	(cudaFree(d_graph_nodes));
	(cudaFree(d_graph_edges));
	(cudaFree(d_graph_mask));
	(cudaFree(d_graph_weights));
	(cudaFree(d_graph_updating_cost));
	(cudaFree(d_cost));
	(cudaFree(d_finished));*/
}
