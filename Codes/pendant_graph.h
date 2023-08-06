#ifndef PENDANT_GRAPH_H
#define PENDANT_GRAPH_H

#include "graph.h"

class PendantGraph {
	public:
	Graph graph;
	vi round;
	vviii vertices_pruned;
	int n_iterations;

	PendantGraph(Graph &g);
	PendantGraph() {}
};
#endif
