#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include <limits.h>
//#include "graph_ear_decompos.h"
#include "ear_graph.h"
#include "bicc.h"
#define BCC_THRESHOLD 40
using namespace std;


void createCSRofG(graph_ear_decompos &ged_obj, int *&g_R, int *&g_F, int *&g_C)
{
	//printf("Graph is\n");
	int num_verts = ged_obj.num_vertices(); 
	int num_edgs = ged_obj.num_edges();
	unsigned int *out_degree_list = ged_obj.get_out_degree_list();
	int *out_array = ged_obj.get_out_vertices();

	g_R = new int[num_verts+1];
	g_C = new int[num_edgs];
	g_F = new int[num_edgs];

	for(int i=0;i < num_verts;i++)
	{
		g_R[i] = out_degree_list[i];

		for(int j = out_degree_list[i]; j < out_degree_list[i+1]; j++)
		{
			g_C[j] = out_array[j];
			g_F[j] = i;
		//	printf("%d %d\n",g_C[j], g_F[j]);
		}	

	}

	//printf("Graph ends\n");
	g_R[num_verts] = num_edgs;

}	

void BiCC::to_eargraph(int g_nodes,int g_edges, vector<vector<pair<pair<int,int>,int > > >&edgelists, map<int,map<int,pair<int,int> > >&bccvertno 
,map<int,int>& artmap, int num_bcc, map<int,int>&bccvrtcnt, map<int,int>& inv_artmap, float *&bc_result, double &g_free_time, int &active_vert_count,
int &free_vert_count)
{

	//Note : artmap has art node and its new value (new value = total num nodes + artpt serial number)
	//Same for inv_art

	int n,m;
	//std::cout<<" In to_eargraph "<<std::endl;
	for(int i=0;i<num_bcc;i++)
	{
		m = 2*edgelists[i].size();	// undirected
		n = bccvrtcnt[i];
		
		if( n <=2 )
			continue;

		//if((( (double)m / (2*g_edges) ) * 100) < 40)
		//	continue;

		std::cout<<"BCC number #: "<<i<<" nodes "<<n<<" edges "<<m<<" Total BCC's are "<<num_bcc<<std::endl;
		
		vector<ear_premble*>ear_premble_vec;
		vector<bool> mark_eargraph_nodes(n,false);
		vector<bool> mark_free_nodes(n,false);
		int *srcs = new int[m];
		int *dsts = new int[m];
		int *wgt  = new int[m];
		int *R, *C, *F;
		int *g_R, *g_C, *g_F;
		int ear_m, ear_n;
		int max_earlevel=0, index = 0;
		ear_premble **ear_info;
		ear_premble **self_ear;
		float *bc_gpu_free = new float[n];

		for(int j = 0 ; j < edgelists[i].size(); j++ )
		{
			int u,v,w;
			srcs[index] = edgelists[i][j].first.first;
			dsts[index] = edgelists[i][j].first.second;
			wgt[index]  = edgelists[i][j].second;
			assert(srcs[index]!=dsts[index]);

		//	cout<<srcs[index]<<" "<<dsts[index]<<" "<<1<<endl;
			index++;
			
			dsts[index] = edgelists[i][j].first.first;
			srcs[index] = edgelists[i][j].first.second;
			wgt[index]  = edgelists[i][j].second;
			index++;
		}
		//exit(0);
		graph_ear_decompos ged_obj;
		ged_obj.init(n,m, srcs, dsts,wgt);   //A001

		//continue;
		bcc_bfs(ged_obj,ear_premble_vec,mark_eargraph_nodes, mark_free_nodes);

		delete [] srcs;
		delete [] dsts;
		delete [] wgt;

		createCSRofG(ged_obj,g_R, g_F, g_C);	
		ged_obj.clear();

		int min_degree_node = INT_MAX;

		vector<int> req_nodes;
	
		for(int k=0;k<n;k++)
		{
			bc_gpu_free[k] = 0;
			if(mark_free_nodes[k]==false)
			{
				req_nodes.push_back(k);
			}
		}	
		//cout<<" Free nodes are "<<req_nodes.size()<<endl;
		if(ear_premble_vec.size() > 0)
			parse_active_nodes(ear_premble_vec, n, R, C, F, ear_m, min_degree_node);

		std::vector<int> parents(n);
		std::vector<int> levels(n,-1);
		std::vector<int> disc_time(n,0);

		if(ear_premble_vec.size() > 0)
			dobfs(R, F, C, n, ear_m, parents,levels,disc_time, min_degree_node);

		for(int k=0;k<levels.size();k++)
		{
			if(max_earlevel<levels[k])
				max_earlevel = levels[k];
		}

		std::vector<int> level_offset(max_earlevel+2,0);
		std::vector<int> level_order(n);
	
		if(ear_premble_vec.size() > 0)
			arrange_level_order(level_offset,level_order,levels,n,max_earlevel,disc_time);
	
		if(ear_premble_vec.size() > 0)
			read_ear_info(ear_premble_vec, R, C, n, ear_m/2, ear_info, self_ear, disc_time, active_vert_count);
	
		disc_time.clear();

		vector<int>artnodes(n,0);
		map<int,pair<int,int> >::iterator ite;
		
		for(ite=bccvertno[i].begin(); ite!=bccvertno[i].end(); ite++)
		{
			if(ite->second.second !=-1)
			{
				artnodes[ite->first] = ite->second.second;
			}	
		}	

		free_vert_count += req_nodes.size();
		if(req_nodes.size()!=0)
		{			
			begin_gpu(g_R, g_C, g_F, n, m, ear_info, self_ear, req_nodes, artnodes, bc_gpu_free, g_free_time);
			for(ite=bccvertno[i].begin(); ite!=bccvertno[i].end(); ite++)
			{
				bc_result[ite->second.first] += bc_gpu_free[ite->first];
			}
		}

		if(ear_premble_vec.size() > 0)
		{
			begin_gpu(g_R, g_C, g_F, n, m, ear_info, parents, levels, level_order, level_offset, max_earlevel, artnodes, bc_gpu_free, g_free_time, R, C);

			for(ite=bccvertno[i].begin(); ite!=bccvertno[i].end(); ite++)
			{
				bc_result[ite->second.first] += bc_gpu_free[ite->first];
			}
		}

		
		if(ear_premble_vec.size() > 0)
		{
			delete [] R;
			delete [] C;
			delete [] F;	
		}	
		
		delete [] g_R;
		delete [] g_C;
		delete [] g_F;
		delete [] bc_gpu_free;
		parents.clear();
		levels.clear();
		level_order.clear();
		level_offset.clear();		
		req_nodes.clear();

	}	
    return ;
}