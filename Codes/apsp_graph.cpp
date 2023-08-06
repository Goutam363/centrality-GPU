/**
 *
 * Parallel APSP implementation using OpenMP
 * Multiple dijkstra's from each possible source
 * Scope : Undirected graphs with positve weights 
 *
 * Input : EdgeList OR adjacency list
 * Output File will contain APSP distance matrix
 **/


#include<iostream>
#include<omp.h>
#include<vector>
#include<queue>
#include<cstdio>
#include<cstdlib>
#include<stack>
#include<cstring>
//#include<pair>
#include "apsp_graph.h"

using namespace std;
vvi APSP_Graph::dist_matrix;
vector<double> APSP_Graph::betweenness;
//vi BiCC::reach;

APSP_Graph::APSP_Graph (int n, viii &edges) {
	graph = Graph(n, edges);
	arr = vi(graph.n_vertices,INT_MAX);
	//	dist_matrix = vvi(n, vi(n, NOT_REACHABLE));
}

APSP_Graph::APSP_Graph (vvii &adj) {
	graph = Graph(adj);
	arr = vi(graph.n_vertices,INT_MAX);
	//	dist_matrix = vvi(adj.size(), vi(adj.size(), NOT_REACHABLE));
}

APSP_Graph::APSP_Graph (Graph &g) {
	//graph = g;
	//arr = vi(graph.n_vertices,INT_MAX);
//	cout<<"graph.n_vertices= "<<g.n_vertices<<" NOT_REACHABLE= "<<NOT_REACHABLE<<endl;	
		dist_matrix = vvi(g.n_vertices, vi(g.n_vertices, NOT_REACHABLE));
		betweenness= vector<double>(g.n_vertices,0);
//		cout<<"Size of dist_matrix[0] = "<<dist_matrix[0].size()<<endl;

	
}

APSP_Graph::APSP_Graph () {
}

void APSP_Graph::parallel_apsp_multiple_dijkstras(int bccnum) {

//int *sigma=new int[graph.n_vertices];
//double *betweenness=new double[graph.n_vertices];
//memset(betweenness,0,sizeof(int)*graph.n_vertices);
//vector<vi> P (graph.n_vertices);		
//stack<int>bet_stack;		

BiCC bicc;
//cout<<" BCC= "<<bccnum<<" vertices= "<<graph.n_vertices<<endl;

//int no_of_nodes = graph.n_vertices;
//int	edge_list_size = graph.n_edges * 2;

int *c_graph_nodes = (int*) malloc(sizeof(int)*graph.n_vertices);
int *c_graph_edges = (int*) malloc(sizeof(int)*graph.n_edges * 2);
//bool *parent_bool = (bool*) malloc(sizeof(bool)*graph.n_edges * 2);
//memset(parent_bool,0,sizeof(bool)*graph.n_edges*2);

int edge_noo = 0;
	
	for ( int i = 0 ; i < graph.n_vertices; i++ ) {
		c_graph_nodes[i] = edge_noo;
		for ( int j = 0; j < graph.adj_list[i].size(); j++ ) {
			c_graph_edges[edge_noo] = graph.adj_list[i][j].first;
		//	h_graph_weights[edge_noo] = graph.adj_list[i][j].second;
		//	cout<<"u = "<<i<<" v= "<<graph.adj_list[i][j].first<<endl;	
			edge_noo++;
		}
	}
//cout<<" nodes= "<<graph.n_vertices<<" edges= "<<graph.n_edges*2<<endl;
//omp_set_num_threads(1);
#pragma omp parallel for 
	for (int i = 0 ; i < graph.n_vertices; i++ ) {

		bool *parent_bool = (bool*) malloc(sizeof(bool)*graph.n_edges * 2);
		memset(parent_bool,0,sizeof(bool)*graph.n_edges*2);
		vi arr1 = vi(graph.n_vertices,NOT_REACHABLE);
		vi sigma= vi(graph.n_vertices,0);
		//cout << i << " : " << omp_get_thread_num() << endl;
		priority_queue<biii, bviii, greater<biii> > q;
		arr1[i] = 0;
		q.push(biii(0,ii(i,i)));  	//q<(distance,(pred,vertex))>
		//cout<<" "<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<" "<<arr[3]<<" "<<arr[4]<<endl;

		//vector<vi> P (graph.n_vertices);		
		stack<int>bet_stack;		
		
		//memset(sigma,0,sizeof(int)*graph.n_vertices);
		//memset(betweenness,0,sizeof(int)*graph.n_vertices);
	
		//cout<<" i is "<<i<<" "<<start<<" "<<num<<" "<<graph.n_vertices<<endl;
		sigma[i]=1;
		while(!q.empty()){
			biii top1=q.top();
			ii top2 = top1.second; 
			q.pop();
			int u = top2.second, d = top1.first;
			int pred=top2.first;			
			bet_stack.push(u);			
			sigma[u]+=sigma[pred];
		//	cout << "sigma "<< u << " " << sigma[u] << endl;
			if ( d <= arr1[u] ) {
				for ( int j = 0 ; j < graph.adj_list[u].size(); j++ ){
					int v = graph.adj_list[u][j].first;
					int w = graph.adj_list[u][j].second;
		//			cout<<" j "<<j<<" u "<<u<<"v "<<v<<" arr[v] "<<arr[v]<<" d+w "<<(d+w)<<endl;
					if(arr1[v] > d + w ) {
			//			cout<<" Test1 u= "<<u <<" v= "<<v <<" size= "<<graph.adj_list[u].size()<<endl;						
						arr1[v] = d + w;
						//dist_matrix[i][v] = d + w;
		//				dist_matrix[ bicc.bccvertno[bccnum][i].first][ bicc.bccvertno[bccnum][v].first ] = d + w;
						q.push(biii(arr1[v],ii(u,v)));
		//				cout<<" Test2"<<endl;						
					//	P[v].push_back(u);
						int ll;
						if(v==graph.n_vertices-1)
							ll=graph.n_edges*2;
						else
							ll=c_graph_nodes[v+1];		
						
						for(int bb=c_graph_nodes[v];bb<ll;bb++)
						{
							if(c_graph_edges[bb]==u)
							{		
								parent_bool[bb]=1;
			//					cout<<" Parent of "<<v<<" is "<<u<<endl;
							}	
						}
						//cout<<" Test3"<<endl;
		//				cout << "v "<<v<<" P[V] "<<P[v].back() <<endl;
						
					}
					else if(arr1[v] == (d + w)){
						//cout<<" Test4"<<endl;							
						sigma[v]+=sigma[u];
					//	P[v].push_back(u);
						int ll;
						if(v==graph.n_vertices-1)
							ll=graph.n_edges*2;
						else
							ll=c_graph_nodes[v+1];		
						
						for(int bb=c_graph_nodes[v];bb<ll;bb++)
						{
							if(c_graph_edges[bb]==u)
							{	
								parent_bool[bb]=1;
			//					cout<<" Parent of "<<v<<" is "<<u<<endl;
							}	
						}
		//				cout << "sigma "<< v << " " << sigma[v] << endl;			
		//				cout << "v "<<v<<" P[V] "<<P[v].back() <<endl;
					}
				
				
				}
			}
		}
	
	/*	for(int i=0;i<graph.n_vertices;i++)
		{
			cout<<"sigma["<<i<<"] "<<sigma[i]<<endl;		
		}
		for(int i=0;i<graph.n_edges*2;i++)
		{
			cout<<"parent_bool["<<i<<"] "<<parent_bool[i]<<endl;		
		}*/
		/*while (!bet_stack.empty())
  		{
		     	cout << " "<< bet_stack.top()<<endl;;
		        bet_stack.pop();
  		}*/
		/*vector<int>::iterator it;
		cout<<"parent = "<<endl;
		for(int i=0;i<graph.n_vertices;i++)
		{
			vi temp=P[i];
			for(it=temp.begin();it!=temp.end();it++)	
				cout<<" "<<*it<<endl;
			
		}*/
	 	accumulate_basic(bet_stack,parent_bool,sigma,i,bccnum,c_graph_nodes,c_graph_edges);
		//P.clear();
	//	break;
		
	}
	/*for(int i=0;i<graph.n_vertices;i++)	
		cout<<" Bet "<<i<<" "<<betweenness[i]<<endl;*/

}

void APSP_Graph::accumulate_basic(stack<int> bet_stack,bool *parent_bool,vi sigma,int s,int bccnum,int *c_graph_nodes,int *c_graph_edges)
{
	BiCC bicc;
 	vector<double> delta (graph.n_vertices);
	delta.assign(graph.n_vertices,0);
	//betweenness[s]+=bet_stack.size()-1;
	int w,v,reach_val=1;
	double coeff;
    //while S:
    while(!bet_stack.empty())
	{	
        	w=bet_stack.top();
		//cout<<"w= "<<w<<endl;
		int artpt=bicc.bccvertno[bccnum][w].first;
		/*if(bicc.artmap.find(artpt)!=bicc.artmap.end())			//To check w is art point or not
			reach_val=bicc.reach[artpt];
		else
			reach_val=1;*/
	//	cout<<" w= "<<w<<" old w= "<<artpt<<endl;	
		if(bicc.bccvertno[bccnum][w].second==-1)
			reach_val=0;
		else
		{
			reach_val=bicc.bccvertno[bccnum][w].second;
			//cout<<"w= "<<w<<" "<<reach_val<<endl;	
			//delta[w]=delta[w]+(2*reach_val);	
		}
		bet_stack.pop();
	        coeff=(1.0+delta[w])/sigma[w];
	        //for v in P[w]:
		//while(!P[w].empty())
		{

			int ll;
			if(w==graph.n_vertices-1)
				ll=graph.n_edges*2;
			else
				ll=c_graph_nodes[w+1];		
				
			for(int bb=c_graph_nodes[w];bb<ll;bb++)
			{
				if(parent_bool[bb]==1)
				{	
					v=c_graph_edges[bb];
				//	cout<<" v= "<<v<<endl; 
					delta[v] += (sigma[v]*coeff);
					delta[v] +=(sigma[v]*coeff*reach_val*2);
		//			cout<<" delta["<<v<<"] = "<<delta[v]<<endl;
		
				}	
			}

			//v=P[w].back();	      
			//delta[v] += (sigma[v]*coeff);
			//delta[v] +=(sigma[v]*coeff*reach_val*2);
			//P[w].pop_back();
		}
//		cout<<endl;        	
		if(w != s)
		{
	            betweenness[artpt] += delta[w];
		//		cout<<" bet= "<<w<<" "<<betweenness[artpt]<<endl;	
		}
		//else
		//	cout<<" bet= "<<w<<" "<<betweenness[artpt]<<endl;		
	}		
/*	cout<<endl;
   	for(int i=0;i<graph.n_vertices;i++)	
		cout<<" delta "<<i<<" "<<delta[i]<<endl;*/
}


/*void APSP_Graph::parallel_apsp_multiple_dijkstras_for_set (int start, int num, int bccnum) {
//#pragma omp parallel for 
//	for ( int i = start ; i < min(start + num, graph.n_vertices); i++ ) 
BiCC bicc;
//cout<<" BCC= "<<bccnum<<" vertices= "<<graph.n_vertices<<endl;
//omp_set_num_threads(2);
//#pragma omp parallel for 
	for (int i = start ; i < min(start + num, graph.n_vertices); i++ ) {
		vi arr1 = vi(graph.n_vertices,NOT_REACHABLE);
		vi sigma= vi(graph.n_vertices,0);
		//cout << i << " : " << omp_get_thread_num() << endl;
		priority_queue<biii, bviii, greater<biii> > q;
		arr1[i] = 0;
		q.push(biii(0,ii(i,i)));  	//q<(distance,(pred,vertex))>
		//cout<<" "<<arr[0]<<" "<<arr[1]<<" "<<arr[2]<<" "<<arr[3]<<" "<<arr[4]<<endl;

		//vector<vi> P (graph.n_vertices);		
		stack<int>bet_stack;		
		
		sigma[i]=1;
		while(!q.empty()){
			biii top1=q.top();
			ii top2 = top1.second; 
			q.pop();
			int u = top2.second, d = top1.first;
			int pred=top2.first;			
			bet_stack.push(u);			
			sigma[u]+=sigma[pred];
		//	cout << "sigma "<< u << " " << sigma[u] << endl;
			if ( d <= arr1[u] ) {
				for ( int j = 0 ; j < graph.adj_list[u].size(); j++ ){
					int v = graph.adj_list[u][j].first;
					int w = graph.adj_list[u][j].second;
		//			cout<<" j "<<j<<" u "<<u<<"v "<<v<<" arr[v] "<<arr[v]<<" d+w "<<(d+w)<<endl;
					if(arr1[v] > d + w ) {
		//				cout<<" Test1"<<endl;						
						arr1[v] = d + w;
						//dist_matrix[i][v] = d + w;
		//				dist_matrix[ bicc.bccvertno[bccnum][i].first][ bicc.bccvertno[bccnum][v].first ] = d + w;
						q.push(biii(arr1[v],ii(u,v)));
		//				cout<<" Test2"<<endl;						
		//				P[v].push_back(u);					
						int ll;
                        if(v==graph.n_vertices-1)
                                ll=graph.n_edges*2;
                        else
                                ll=c_graph_nodes[v+1];
                        for(int bb=c_graph_nodes[v];bb<ll;bb++)
                        {
                                if(c_graph_edges[bb]==u)
                                {
                                        parent_bool[bb]=1;
                                       // cout<<" Parent of "<<v<<" is "<<u<<endl;
                                }
                        }

						//cout<<" Test3"<<endl;
		//				cout << "v "<<v<<" P[V] "<<P[v].back() <<endl;
						
					}
					else if(arr1[v] == (d + w)){
						//cout<<" Test4"<<endl;							
						sigma[v]+=sigma[u];
			//			P[v].push_back(u);
						int ll;
                        if(v==graph.n_vertices-1)
                           ll=graph.n_edges*2;
                        else
                           ll=c_graph_nodes[v+1];
                        for(int bb=c_graph_nodes[v];bb<ll;bb++)
                        {
                                if(c_graph_edges[bb]==u)
                                {
                                        parent_bool[bb]=1;
                                       // cout<<" Parent of "<<v<<" is "<<u<<endl;
                                }
                        }
			}
		}
	 	//accumulate_basic(bet_stack,P,sigma,i,bccnum);
		P.clear();
		break;
		
	}


}*/
void APSP_Graph::print_dist_matrix_to_file(char *file) {
	FILE *fp = fopen(file, "w");
	if ( !fp ) 
		cout << "Error in opening file for writing\n", exit(3);

	for(int i= 0 ; i < graph.n_vertices; i++ ){
		for (int j =0; j < graph.n_vertices; j++ )
			fprintf(fp, "0 ");
		fprintf(fp, "\n");
	}
}

void APSP_Graph::print_dist_matrix() {
	for(int i= 0 ; i < graph.n_vertices; i++ ){
		for (int j =0; j < graph.n_vertices; j++ )
			continue;
		//printf("%d ", dist_matrix[i][j]);
		printf("\n");
	}
}

int APSP_Graph::get_max_distance() {
	vi dist_v(graph.n_vertices,0);

#pragma omp parallel for
	for(int i= 0 ; i < graph.n_vertices; i++ ) {
		int max_dist = 0;
		for (int j =0; j < graph.n_vertices; j++ )
			continue;
		//		max_dist = max(max_dist, dist_matrix[i][j]);
		dist_v[i] = max_dist;
	}

	int max_dist = 0;
	for (int i=0; i<graph.n_vertices; i++) 
		max_dist = max(max_dist, dist_v[i]);
	return max_dist;
}
/*
   void APSP_Graph::parallel_apsp_floyd_warshal_compare () {
   vvi dist = vvi ( graph.n_vertices, vi (graph.n_vertices, NOT_REACHABLE));
   for ( int i = 0 ; i < graph.n_vertices; i++ ) {
   dist[i][i] = 0;
   for ( int j = 0 ; j < graph.adj_list[i].size(); j ++ ){
   int v = graph.adj_list[i][j].first;
   int w = graph.adj_list[i][j].second;
//dist[i][v] = w;
//dist[v][i] = w;
}
}

for ( int k = 0 ; k < graph.n_vertices; k ++ ) {
#pragma omp parallel for
for ( int i = 0 ; i < graph.n_vertices; i++ ) {
for ( int j = 0 ; j < graph.n_vertices; j++ ) {
if ( dist[i][k] != NOT_REACHABLE && dist[k][j] != NOT_REACHABLE) 
dist[i][j] = min( dist[i][j], dist[i][k] + dist[k][j]);
}
}
}

bool same = true;
for ( int i = 0 ; same && i < graph.n_vertices; i++ )
for ( int j = 0 ; same && j < graph.n_vertices; j++ )
if ( dist[i][j] != dist_matrix[i][j] )
same = false;

for(int i= 0 ; i < n_vertices; i++ ){
for (int j =0; j < n_vertices; j++ )
cout << dist[i][j] << " ";
cout<<endl;
}

cout << "Same : " << (same ? "true" : "false") << endl;
}
*/
