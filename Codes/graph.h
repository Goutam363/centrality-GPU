 /** Using Adjancency List :
 * adj_list[i][j] contains j'th neighbour and weight of the edge to that neighbour
 *	Init  	: 	Graph graph = Graph ( filename);
 *			OR
 *			Graph graph = Graph ( int n_vertices, vector<pair<pair<int,int>,int>> edgelist );
 *			OR
 *			Graph graph = Graph ( vector<vector<pair<int,int>>> adjacency_list );
 *
 **/

#ifndef GRAPH_H_
#define GRAPH_H_

#include<vector>

using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef pair< ii, int> iii;
typedef pair< ii, ii> iiii;
typedef vector<ii> vii;
typedef vector<iii> viii;
typedef vector<viii> vviii;
typedef vector<iiii> viiii;
typedef vector<viiii> vviiii;
typedef vector<vii> vvii;
typedef pair<int,ii> biii;
typedef vector<biii> bviii;
#define mp(x,y) make_pair(x,y)
#define pb(x) push_back(x)

class Graph {
	public :
		int n_vertices;
		int n_edges;
		vvii adj_list;

		Graph();
		Graph(int n, viii &edges);
		Graph(vvii &adj);
		Graph(char *file);
};
#endif 
