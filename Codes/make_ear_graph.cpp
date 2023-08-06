#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>
//#include "graph_ear_decompos.h"
#include "ear_graph.h"

using namespace std;

void dfs(vector<bool>& visited, map<int,vector<int> >& adj, int start, int& last_v)
{

//	cout<<" stat "<<start<<endl;

	visited[start] = true;
	
	for(int i=0;i<adj[start].size();i++)
	{
		if(visited[adj[start][i]]==false)
		{
			dfs(visited, adj, adj[start][i], last_v);
			if(last_v==-1)
				last_v = adj[start][i];
		}	
	}	
}


void connect_graphs(map<int,vector<int> >& adj, int n, vector<bool>& visited, int &edge_counter)
{
	
	int start=0,block=0;
	for(int i=0;i<n;i++)
	{
		if(adj.find(i)==adj.end())
		{
			visited[i]= true;
		}else
		{
			if(block==0)
			{
				start = i;
				block = 1;	
			}	
			
		}	
	}		

	int last_v = -1;
	dfs(visited,adj,start, last_v);
	
	for(int i=0;i<n;i++)
	{
		if(visited[i]==false)
		{
			adj[last_v].push_back(i);
			adj[i].push_back(last_v);
			edge_counter += 2;
			last_v = -1;
			dfs(visited,adj,i, last_v);
			
		}	
	}	

}

void create_CSR(map<int,vector<int> >& adj, int*&R, int *&C, int *&F, int m, int n)
{

	R = new int[n+1];
	C = new int[m];
	F = new int[m];

	for(int i=0;i<n;i++)
	{
		R[i] = 0;
		if(adj.find(i)!=adj.end())
		{
			R[i] = adj[i].size();
		}	
	}

	R[n] = m;

	int t = R[0];				// 0 1 2 3 4
	R[0] = 0;					// 5 3 2 4 2
	for(int i=0;i<n-1;i++)		// 0 5 8 10 14 16
	{
		int d = R[i+1];
		R[i+1] = R[i] + t;
		t = d;
	}

	assert((t+R[n-1])==R[n]);

	int k=0;
	for(int i=0;i<n;i++)
	{

		if(adj.find(i)!=adj.end())
		{
			k=0 ;
			for(int j=R[i]; j<R[i+1]; j++)
			{
				F[j] = i;
				C[j] = adj[i][k];
				k++;

			}
			//cout<<" "	
			assert(adj[i].size()==k);
		}	
	}	
}

void parse_active_nodes(vector<ear_premble*>&ear_premble_vec, int n, int *&R, int *&C, int*&F, int &ear_m, int &min_degree_node)
{

	map<int,vector<int> > adj;
	int edge_counter = 0, index=0;
	int degree_2_vert_cnt = 0;
	set<int> active_vertices;
	for(int kk=0;kk<ear_premble_vec.size();kk++)
	{

		int src = ear_premble_vec[kk]->lnode;
		int dst = ear_premble_vec[kk]->rnode;

		if(ear_premble_vec[kk]->cnt > 0)
			degree_2_vert_cnt += ear_premble_vec[kk]->cnt;	
		
		//cout<<"src "<<src<<" dst "<<dst<<endl; 
		active_vertices.insert(src);
		active_vertices.insert(dst);

		if(src==dst)
		{
			//self_loop.insert(src);
			index++;
			continue;
		}	
		//std::cout<<src<<" "<<dst<<" "<<std::endl;

		int flag = 0;
				
		if(adj.find(src)==adj.end())
		{
			std::vector<int> tmp;
			tmp.push_back(dst);
			adj[src] = tmp;
			++edge_counter;
		}
		else
		{	
			for(unsigned int i=0;i<adj[src].size();i++)
			{
				if(adj[src][i]==dst)
				{
					flag = 1;
				}
			}
			if(flag==0)
			{	
				adj[src].push_back(dst);
				++edge_counter;
			}
		}	

			

		flag = 0;
		if(adj.find(dst)==adj.end())
		{
			std::vector<int> tmp;
			tmp.push_back(src);
			adj[dst] = tmp;
			++edge_counter;
		}
		else
		{
			for(unsigned int i=0;i<adj[dst].size();i++)
			{
				if(adj[dst][i]==src)
				{
					flag = 1;
				}
			}	
			if(flag==0)
			{	
				adj[dst].push_back(src);
				++edge_counter;
			}
		}

		if(adj[src].size() > 0 && adj[src].size() < min_degree_node)	
			min_degree_node = src;
		if(adj[dst].size() > 0 && adj[dst].size() < min_degree_node)
			min_degree_node = dst;

	}

/*	cout<<" vertices in ear graph "<<adj.size()<<" edges are "<<edge_counter<<endl;	
	cout<<" dependent vertices count is "<<active_vertices.size()<<endl;;
	cout<<" degree 2 vertices count is "<<degree_2_vert_cnt<<endl;;*/

	vector<bool>visited(n,false);
	connect_graphs(adj,n,visited,edge_counter);
	visited.clear();
	create_CSR(adj, R, C, F, edge_counter, n);
	ear_m = edge_counter;
	adj.clear();

	assert(index==0);
	active_vertices.clear();
}
