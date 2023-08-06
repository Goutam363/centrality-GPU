/**
 * @author : Ashutosh Kumar ( ashumac@gmail.com)
 *
 * Implementation of Graph Class.
 **/

#include"graph.h"
#include<cstdlib>
#include<cstdio>
#include<iostream>

Graph::Graph() {
	n_vertices = n_edges = 0;
}

Graph::Graph(int n, viii &edges) {
	n_vertices = n;
	n_edges = edges.size();
	adj_list = vvii(n_vertices);
	for(int i = 0 ; i < n_edges; i++) {
		int u = edges[i].first.first;
		int v = edges[i].first.second;
		int w = edges[i].second;
		adj_list[u].pb( mp(v,w) );
		adj_list[v].pb( mp(u,w) );
	}
}

Graph::Graph(vvii &adj) {
	adj_list = adj;
	n_vertices = adj.size();
	n_edges = 0;
	for(int i=0; i < n_vertices; i++ )
		n_edges += adj[i].size();
	n_edges /= 2;
}

Graph::Graph(char *file) {
	FILE *fp = fopen(file, "r");
	if (!fp) {
		cout << "Input file cannot be opened\n";
		exit(2);
	}
	fscanf(fp, "%d %d", &n_vertices, &n_edges);
	
	adj_list = vvii(n_vertices);
	vi degree_2;//= vi(n_vertices);
	degree_2.assign(n_vertices,0);

	for(int i = 0 ; i < n_edges ; i++ ){
		int u, v, w;
		fscanf(fp, "%d %d %d", &u, &v, &w);
		//u--;
		//v--;
		adj_list[u].pb( mp(v,w) );
		degree_2[u]++;	
		adj_list[v].pb( mp(u,w) );
		degree_2[v]++;	
	}

	int fr_cnt=0;
	int two_cnt=0;
	for(int i=0;i<n_vertices;i++)
	{
		if(degree_2[i]==4)
			fr_cnt++;
		if(degree_2[i]==2)
			two_cnt++;
	}
//	cout<<" 2 n 4 degree count "<<two_cnt<<" "<<fr_cnt<<endl;
//	exit(0);
}
