/**
 * @author : Ashutosh Kumar (ashumac@gmail.com)
 * 
 * Class strucutre for APSP based on Biconnected Components.
 *	Partitioning the graph based on the BiCC, then running our hybrid algorithm using both CPU and GPU.
	
 *	Init 	: 	Modified_APSP mod_apsp;
 *	Usage 	: 	mod_apsp.run_modified_apsp(graph, USE_CPU, USE_GPU, REMOVE_PENDANTS);
 *	Output	: 	Write to file / Find Maximum Distance / Print
 *
 **/

#ifndef MODIFIED_APSP_H_
#define MODIFIED_APSP_H_

#include "pendant_graph.h"
//#include "apsp_graph.h"
#include <map>
#include <utility>


class Modified_APSP {
	public : 
	//	vector<APSP_Graph> apsp_graphs;
		vii new_id;
		//vector<Graph> graphs;
		vviiii shortest_ap;
		//APSP_Graph condensed_graph;
		int n_vertices;
		int n_bcc;
		vii bcc_sizes;
		int gpu_graphs;
		bool use_cpu;
		bool use_gpu;
		bool remove_pendants;
		PendantGraph p_graph;
		//float *bc_result;

		Modified_APSP () {}
		void run_dummy_memcpy();
		
		void run_modified_apsp(Graph &graph, bool cpu, bool gpu, bool pendants, float *&bc_result, double &g_free_time, int &active_vert_count,
		int &free_vert_count);

		void print_dist_matrix_to_file (char *file);
		int get_max_distance();
		void *cpu_worker();
		void *gpu_worker();
		
};
#endif
