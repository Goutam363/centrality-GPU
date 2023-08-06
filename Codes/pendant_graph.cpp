#include"pendant_graph.h"
#include<climits>
#include<iostream>
#include<map>

using namespace std;

#define ANATH INT_MAX
PendantGraph::PendantGraph(Graph &g) {
	round = vi(g.n_vertices,-1);
	n_iterations = 0;
	for ( int i = 0; ; i++ ) {
		viii iteration;
		for ( int j = 0; j < g.n_vertices; j++ ) {
			if ( round[j] != -1 )
				continue;
			int nbrs = 0;
			ii par = mp(ANATH,ANATH);
			for ( int k = 0 ; k < g.adj_list[j].size(); k++ ) {
				if ( round[g.adj_list[j][k].first] == -1 ){
					par = g.adj_list[j][k];
					nbrs ++;
				}
			}

			if ( nbrs <= 1 )
				iteration.pb(mp(par,j));
		}
		for ( int j = 0 ; j < iteration.size(); j++ ) {
			round[iteration[j].second] = i;
		}
		n_iterations++;
		//cout << "Vertices removed in iteration " << n_iterations << " : " << iteration.size() << endl;
		if ( (iteration.size() ==0 or i>=4 ) ) /* * 100.0)/g.n_vertices  < 0.05 )*/
			break;
		vertices_pruned.pb(iteration);
	}

	vi new_id(g.n_vertices, -1);
	int cnt = 0;
	for (int i = 0; i < g.n_vertices; i++ ) {
		if ( round[i] == -1 ) {
			new_id[i] = cnt;
			cnt++;
		}
	}
	vvii adj;
	for (int i = 0; i < g.n_vertices; i++ ) {
		if ( round[i] == -1 ) {
			vii list;
			for ( int j = 0 ; j < g.adj_list[i].size(); j++ )
				if ( round[g.adj_list[i][j].first] == -1 )
					list.pb(ii(new_id[g.adj_list[i][j].first],g.adj_list[i][j].second));
			adj.pb(list);
		}
	}
	graph = Graph(adj);
}

