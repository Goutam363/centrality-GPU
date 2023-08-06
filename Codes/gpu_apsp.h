/**
 * (Based on work of Pawan Harish)
 *
 * Class structure of APSP code on GPU. 
 *	Perform parallel APSP on GPU usign multiple SSSP. 
 *
 *	Init 	: 	GPU_APSP gpu_apsp;
 *			gpu_apsp.set_graph(graph);
 *	Usage	: 	gpu_apsp.run_apsp(apsp_graph);
 *	Free	:	gpu_apsp.free_graph();
 **/

#ifndef GPU_APSP_H
#define GPU_APSP_H

#include "graph.h"
#include "apsp_graph.h"

class GPU_APSP {
	int no_of_nodes;
	int edge_list_size;
	int num_of_blocks;
	int num_of_threads_per_block;
	bool finished;
	
	int *h_graph_nodes;
	int *h_graph_edges;
	int *h_cost;
	bool *h_graph_mask;
	int *h_graph_updating_cost;
	short int *h_graph_weights;
	int *h_old_vertex;
	int *h_new_vertex;
	int *h_reachval;
	
	double *d_betweenness;
	double *h_betweenness;	
	
	int* d_graph_nodes;
	int* d_graph_edges;
	int* d_cost;
	bool* d_graph_mask;
	int* d_graph_updating_cost;
	short int* d_graph_weights;
	bool *d_finished;
	//int *d_upq,*d_vpq;
	int *dist;	
	int *d_old_vertex;
	int *d_new_vertex;
	int *d_reachval;


	int *h_arr1;
	int *h_sigma;
	int *h_dist;
	int *h_upq,*h_vpq;
	int *h_pq_index;
	int *h_seq;
	int *h_bstack;
	int *h_bstackind;
	int *h_pindex;	
	int *h_t_parent;
	double *h_delta;	
	bool *h_parent_bool;

	int *d_arr1;
	int *d_sigma;
	int *d_dist;
	int *d_upq,*d_vpq;
	int *d_pq_index;
	int *d_seq;
	int *d_bstack;
	int *d_bstackind;
	int *d_pindex;
	int *d_t_parent;
	double *d_delta;		
	bool *d_parent_bool;

	public :
	void run_dummy_memcpy();
	void set_graph(Graph &graph,int bccnum);
	void free_graph();
	void run_apsp(int start, int num,int totalvertices,int bccnum );
	//void run_apsp(APSP_Graph &graph);
	//void run_apsp();
};
#endif
