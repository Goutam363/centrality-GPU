/**

**/

#include<iostream>
#include"bicc.h"
#include<algorithm>
#include<cstdio>
#define INVALID -1


using namespace std;
int children=0;
//int sumbccvertcnt=0;
BiCC::BiCC (char *file) {
	graph = Graph(file);
	n_vertices = graph.n_vertices;
	last_dfs_number = 0;
	bcc_cnt=-1;
	prv_u=-1;	
	dfs_number = vi(n_vertices, INVALID);
	high_water = vi(n_vertices, INVALID);
}

BiCC::BiCC (Graph &g) {
	graph = g;
	n_vertices = graph.n_vertices;
	last_dfs_number = 0;
	bcc_cnt=-1;
	prv_u=-1;	
	dfs_number = vi(n_vertices, INVALID);
	high_water = vi(n_vertices, INVALID);
	block_index.assign(graph.n_vertices,mp(0,0));
	//reach= vi(graph.n_vertices,0);
}
BiCC::BiCC()
{
}
vviii BiCC::run_bicc() {
	//cout<<" The total vertices are "<<n_vertices<<endl;
	for ( int i=0 ; i < graph.n_vertices; i++){
		if (dfs_number[i] == INVALID )
			bicc_tarjan(i, INVALID);
			}
	/*cout<<" The total vertices are "<<n_vertices<<endl;
	for( int i=0 ; i < graph.n_vertices ; i++ )
		{
		cout<<"charu'i is "<<i<<" low is "<<high_water[i]<<endl;  
		}*/	
	return subgraphs_edgelist;
}

void BiCC::print_bicc_hist_to_file(char *file) {
	n_vertices = graph.n_vertices;
	int n_bcc = subgraphs_edgelist.size();
	FILE *fp = fopen(file, "w");

	if ( n_bcc == 0 ) 
		return;

	vii bcc_sizes;
	for (int i = 0; i < n_bcc; i++ )
		bcc_sizes.pb(mp(subgraphs_edgelist[i].size(), i));
	sort(bcc_sizes.rbegin(), bcc_sizes.rend());

	int prev = bcc_sizes[0].first;
	int bcc_cnt = 1;
	for ( int i=1; i <= n_bcc; i++ ) {
		if ( (i < n_bcc) and (bcc_sizes[i].first == prev )) 
			bcc_cnt++;
		else {
			fprintf(fp, "%d %d\n", prev, bcc_cnt);
			bcc_cnt = 1;
			if ( i < n_bcc )
				prev = bcc_sizes[i].first;
		}
	}
	fclose(fp);
}

void BiCC::bicc_tarjan(int u, int parent ) {
	last_dfs_number++;
	int b1=0;
	dfs_number[u] = last_dfs_number;
	high_water[u] = last_dfs_number;
	if(parent==0)
		children++;
	for ( int i= 0; i < graph.adj_list[u].size(); i++ ) {
		int v = graph.adj_list[u][i].first;
		//cout<<"u is ::"<<u<<""<<"v is::"<<v<<endl;
		int w = graph.adj_list[u][i].second;
		if ( dfs_number[v] == INVALID ) {
			edge_stack.push(mp(mp(u,v),w));
			bicc_tarjan(v, u);
			if ( high_water[v] < high_water[u])
				high_water[u] = high_water[v];
			if ( dfs_number[u] <= high_water[v] && edge_stack.size()!=0) {
				
				int a, b, w,key;
				map <int,bool> bccvertcnt;	//To find numberof vertices in bcc;
				viii edge_list;
				++bcc_cnt;
				//int flag=0;		
				int tempstr;
				//cout<<"Test 1 :: "<<u<<" "<<v<<endl;//<<" Addr edge_list "<<edge_list<<endl;		//CH01	
					
				a = edge_stack.top().first.first;
				if(bccvertcnt.find(a)==bccvertcnt.end())
					bccvertcnt.insert(mp(a,0));
				block_index[a]=mp(bcc_cnt,0);
				//if(a==prv_u)
				//	flag=1;
				b = edge_stack.top().first.second;
				if(bccvertcnt.find(b)==bccvertcnt.end())
					bccvertcnt.insert(mp(b,0));	
				block_index[b]=mp(bcc_cnt,0);
				//if(b==prv_u)
				//	flag=1;
				edge_list.pb(edge_stack.top());
				edge_stack.pop();
				while (!(a == u && b == v )) {
					if(edge_stack.size()==0)
						break;
					a = edge_stack.top().first.first;
					if(bccvertcnt.find(a)==bccvertcnt.end())
						bccvertcnt.insert(mp(a,0));
					block_index[a]=mp(bcc_cnt,0);
				//	if(a==prv_u)
				//		flag=1;					
					b = edge_stack.top().first.second;
					if(bccvertcnt.find(b)==bccvertcnt.end())
						bccvertcnt.insert(mp(b,0));
					block_index[b]=mp(bcc_cnt,0);
				//	if(b==prv_u)
				//		flag=1;
					edge_list.pb(edge_stack.top());
					edge_stack.pop();
				}
				subgraphs_edgelist.pb(edge_list);	//Starts from 0
				bccmapping.pb(bcc_cnt);			//Starts from 1..holds mapping of BCC's to integer value
				if((u!=0) || (u==0 && children > 1))
					if(artmap.find(u)==artmap.end())
					{
						vertices_count++;
						artmap.insert(pair<int,int>(u,vertices_count));
						//sumbccvertcnt+=bccvertcnt.size();		
						//reach[u]+=sumbccvertcnt;				//Calculate reachability
						//cout<<" reach["<<u<<"]= "<<reach[u]<<endl;
						inv_art.insert(pair<int,int>(vertices_count,u));
						//cout<<" Art u= "<<u<<" Alias vert= "<<vertices_count<<endl;				
					}
							
				//	prv_u=u;
			} 
		} else if ( dfs_number[v] < dfs_number[u] && v != parent ) {
			edge_stack.push(mp(mp(u,v),w));
			if ( dfs_number[v] < high_water[u]) 
				high_water[u] = dfs_number[v];
		}
	}
}
