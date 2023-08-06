#ifndef HYBRID_APSP_
#define HYBRID_APSP_

#include "apsp_graph.h"
#include "gpu_apsp.h"
#include "pendant_graph.h"

//#include<atomic>

#define CHUNK_SIZE 1000
#define GPU_MIN_SIZE 1000

using namespace std;

class Hybrid_APSP {
	PendantGraph p_graph;
	APSP_Graph apsp_graph;
	GPU_APSP gpu_apsp;
	bool use_gpu;
	bool use_cpu;
	bool remove_pendants;
	int n_vertices;
	int n_chunks;
//	atomic_int gpu_chunk;
//	atomic_int cpu_chunk;
	int gpu_chunk;
	int cpu_chunk;

	public : 
		Hybrid_APSP(Graph &g, bool cpu, bool gpu, bool pendant);
		void run_hybrid_apsp();
		void * cpu_worker();
		void * gpu_worker();
		void print_dist_matrix_to_file(char *);
};
#endif 
