/**
 *
 * Complie : sh compile.sh 
 *
 * RUN : ./apsp INPUT_FILE_NAME OUTPUT_FILE_NAME [NUM_THREADS]
 * Input File : Should have edge list.
 * Output File will contain APSP distance matrix
 **/

#include "apsp_graph.h"
#include "bicc.h"
#include "modified_apsp.h"
#include "pendant_graph.h"
#include "hybrid_apsp.h"
#include <sys/time.h>
#include <time.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <map>
#include <cstring>

#define N_LIMIT 20001

double run_modified_apsp (Graph &graph, bool cpu, bool gpu, bool pendants, char *out_file, double &g_free_time) {

	//int max_dist;
	struct timeval start_time, end_time;
	gettimeofday(&start_time, NULL);

	float *bc_result = new float[graph.n_vertices];
	int active_vert_count = 0;
	int free_vert_count = 0;
	for(int i=0; i < graph.n_vertices;i++)
	{
			bc_result[i] = 0;
	}	

	Modified_APSP mod_apsp = Modified_APSP();
	mod_apsp.run_modified_apsp(graph, cpu, gpu, pendants, bc_result, g_free_time, active_vert_count, free_vert_count);

	double run_time = 2.0;

	/*printf("Betweenness Centrality is ============\n");
	for(int i=0; i < graph.n_vertices;i++)
	{
			printf("%d %f\n",i,bc_result[i]);
	}*/
	printf("active_vert_count %d free_vert_count %d \n",active_vert_count,free_vert_count );		
	return run_time;
}

int main (int argc, char **argv){
	if ( argc < 2 )
		printf("Usage : ./apsp INPUT_FILE_NAME1 INPUT_FILE_NAME2 .. [OUTPUT_FILE_NAME] [NUM_THREADS]\n"), exit(1);

	char *filename_prefix = argv[1];
	
	double g_free_time=0;
	//for ( int i = 1; i < argc; i++ )
	{
		Graph graph = Graph(argv[1]);

		double mod_hybrid_time = run_modified_apsp(graph, true, true, false, filename_prefix, g_free_time);	
		//printf("%s,%.3lf\n", argv[i],mod_hybrid_time);
		
	}

	printf("GPU total time is %f \n",g_free_time);
	return 0 ;
}
