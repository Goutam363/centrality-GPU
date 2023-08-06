/**
 * @author : Ashutosh Kumar ( ashumac@gmail.com )
 *
 * Class Definition of the APSP related functions.
 *Storing the Graph Object as well as N*N distance matrix for APSP.
 *	Init 	:	Similar to Graph.	
 *	Running :	apsp_graph.parallel_apsp_multiple_dijkstras();
 *			OR
 *	 		apsp_graph.parallel_apsp_floyd_warshal_compare();
 *	
 *	Output 	:	Write to file / Find maximum distance / Print.
 **/

#ifndef APSP_GRAPH_H_
#define APSP_GRAPH_H_

#include"graph.h"
#include<climits>
#include"bicc.h"

#define NOT_REACHABLE INT_MAX

class APSP_Graph {
	public :
	Graph graph;
	static vvi dist_matrix;
	static vector<double> betweenness;
	vi arr;
	APSP_Graph (int n, viii &edges);
	APSP_Graph (vvii &adj);
	APSP_Graph (Graph &g);
	APSP_Graph ();
	void print_dist_matrix_to_file(char *file);
	void print_dist_matrix();
	int get_max_distance();
	void parallel_apsp_multiple_dijkstras(int a);
	void parallel_apsp_multiple_dijkstras_for_set(int a, int b,int bccnum);
	void parallel_apsp_floyd_warshal_compare();
	void accumulate_basic(stack<int> bet_stack,bool *parent_bool,vi sigma,int s,int bccnum,int *c_graph_nodes,int *c_graph_edges);
};
#endif
