/**
 * @author : Ashutosh Kumar (ashumac@gmail.com)
 *
 * Implementation of Our Algorithm
 **/

#include "modified_apsp.h"
#include "bicc.h"
#include<cstdio>
#include<omp.h>
#include<algorithm>
#include<iostream>
#include "gpu_apsp.h"
#include<sys/wait.h>
#include<sys/types.h>
#include "utils.h"
#include "hybrid_apsp.h"
#include<cstddef>
#include <time.h>
#include <sys/time.h>
#include<queue>
#include <list>
#include <fstream>
map<int,map<int,pair<int,int> > > BiCC::bccvertno;	//block new vertex old vertex mapping
map<int,map<int,int> > BiCC::rev_bccvertno;		//Block oldartpt to new_art mapping <bccnumber,<relabeldval,oldvertex>>
vector<int> BiCC::blocklevel;		//To store BCC levels
int avglevel1=0;
int avglevel2=0;
int maxlevel=0;
int templevel=0;
int newbcccnt= 0;
int CPU_WORK_UNIT_SIZE =100;
int WORK_UNIT_SIZE=6000;
using namespace std;
int cpu_chunk,gpu_chunk;
int var,cpu_var;
int flag=0;
long int total_vertices,total_bcc;


void BiCC::print_edgelist()
{
	ofstream f1;
	f1.open("bccmap.txt",ios::out);
	f1<<" blocklevel vector size= "<<blocklevel.size()<<endl;
	for(int i=0;i<blocklevel.size();i++)
	{
		f1<<"i= "<<i<<" block_level= "<<blocklevel[i]<<endl;
	}
	map<int,viii>::iterator itr;
	map<int,map<int,int> >::iterator mitr;	
	map<int,int>::iterator mitr2;
	cout<<" bcc map "<<endl;
	map<int,vi>::iterator mitr3;
	vi temp;
	for(mitr3=bccmap.begin();mitr3!=bccmap.end();mitr3++)
	{
		temp=mitr3->second;
		f1<<" Art= "<<mitr3->first<<" "<<"list:  ";
		for(int j=0;j<temp.size();j++)
		{
			f1<<temp[j]<<" , ";
		}
		f1<<endl;
	}	
	f1<<" Inv art is "<<endl;
	//map<int,int>::iterator mitr2;
	for(mitr2=inv_art.begin();mitr2!=inv_art.end();mitr2++)	
	{
		//if(mitr2)
		f1<<"new vertex= "<<mitr2->first<<" old vertex= "<<mitr2->second<<endl;	
	}


}

void BiCC::fillbccparent(int totalvertex)
{
	vector<bool> visited;
	visited.assign(vertices_count+1,false);
	bccparent.assign(vertices_count+1,-1);
	blocklevel.assign(vertices_count+1,-1);
	leafnodes.assign(vertices_count+1,false);
	artdegree.assign(vertices_count+1,0);
	list<int> queue;
	int level=0;
	visited[0]=true;
	queue.pb(0);
	int mark=0;
	int mark_initial=0;
	int flag=0;
	int popped_index=-1;
	vi::iterator i;
	vi tempvector;		
	vi::iterator vitr;
	blocklevel[0]=level;
	int curr_level=0;
	level++;
	ofstream f2;
	f2.open("parent.txt");
	while(!queue.empty())
	{
		int s=queue.front();
		popped_index++;
		queue.pop_front();
		tempvector=(bccmap.find(s)->second);

		for(i=tempvector.begin();i!=tempvector.end();i++)
		{
			//cout<<" "<<*i<<endl;
			if(!visited[*i])
			{
				visited[*i]=true;
				bccparent[*i]=s;
				artdegree[*i]+=1;
				artdegree[s]+=1;
				leafnodes[s]=true;	
				//f2<<*i<<"   "<<s<<endl;	
				blocklevel[*i]=level;
				if(level> maxlevel)
					maxlevel=level;
				if((*i)<totalvertex)
				{
					if(bcclevelmap.find(level)!=bcclevelmap.end())
					{
						//vi tempv=bcclevelmap[level];
						bcclevelmap[level].pb(*i);		
					}
					else
					{
						vi tempv;
						tempv.pb(*i);
						bcclevelmap[level]=tempv;
						
					}
				}
				queue.push_back(*i);
				flag++;
			}
		}
		if(mark_initial = 0)
		{
			mark=flag;
			mark_initial = 1;
			level++;
		}
		if(popped_index == mark)
		{
			mark = flag;
			level++;
		}		
		
		
	}	
}

void BiCC::reachfunction(int totalvertex, float *&bc_result)
{
	//APSP_Graph ag;
	int b,bccnum,artpt,rel_vrt,cnt=0;
	vector<bool> visitedbcc;
	vi reachvec;
	int flag=0;
	visitedbcc.assign(newbcccnt,false);
	reachvec.assign(vertices_count+1,0);

	int templevel=maxlevel;
	while(maxlevel!=0)
	{
		vi tempv=bcclevelmap[maxlevel];
		for(int tempvindex=0;tempvindex<tempv.size();tempvindex++)
		{	
				b=tempv[tempvindex];
				cnt=0;
				cnt=reachvec[b]+bccvrtcnt[b]-1;
				bccnum=b;
				b=bccparent[b];

				if(b>totalvertex)
				{
					int cntx,cnty;
					artpt=inv_art[b];
					rel_vrt=rev_bccvertno[bccnum][artpt];
					if(bccvertno[bccnum][rel_vrt].second==-1)
					{
						bccvertno[bccnum][rel_vrt].second=(totalvertex-cnt-1);
						cntx=bccvertno[bccnum][rel_vrt].second;
					}					
					else
					{
						bccvertno[bccnum][rel_vrt].second=(bccvertno[bccnum][rel_vrt].second-cnt);
					}
					
					bccnum=bccparent[b];
					reachvec[bccnum]+=cnt;
					rel_vrt=rev_bccvertno[bccnum][artpt];
					if(bccvertno[bccnum][rel_vrt].second==-1)
					{			
						bccvertno[bccnum][rel_vrt].second=cnt;
						cnty=bccvertno[bccnum][rel_vrt].second;					
					}					
					else
					{
						bccvertno[bccnum][rel_vrt].second+=cnt;
						cntx=cnt;
						cnty=(totalvertex-bccvertno[bccnum][rel_vrt].second-1);
					}
					
						bc_result[artpt]+=(cntx*cnty*2); //To Match with baders results
					flag=0;
				}
				else
				{	
					cout<<bccnum<<endl;
				}
				
		}	
		maxlevel-=2;
	}
}

//To print final N x N  distance matrix
void BiCC::printoutput(float *&bc_result)
{
	cout<<" BetCent output  "<<endl;
	for(int i=0;i<graph.n_vertices;i++)
	{
		cout<<bc_result[i]<<" ";	
	}
	cout<<endl;
}
void Modified_APSP::run_modified_apsp(Graph &graph, bool cpu, bool gpu, bool pendants, float *&bc_result, double &g_free_time, int &active_vert_count,
int &free_vert_count ) {

	struct timeval start_time, end_time;
	double run_time;
	gettimeofday(&start_time, NULL);

	use_cpu = cpu;
	use_gpu = gpu;
	remove_pendants = pendants;

	if ( remove_pendants ) {
		p_graph = PendantGraph(graph);
		graph = p_graph.graph;
	}
	//cout << "graph size " << graph.n_vertices << " vertices and " << graph.n_edges << " edges." << endl;
	BiCC bicc = BiCC(graph);
	bicc.vertices_count=graph.n_vertices;
	vviii edgelists = bicc.run_bicc();
	
	n_vertices = graph.n_vertices;
	//APSP_Graph ag=APSP_Graph(graph);
	gettimeofday(&end_time, NULL);
	run_time = ((end_time.tv_sec - start_time.tv_sec)*1000 + (end_time.tv_usec - start_time.tv_usec)/1000.0)/1000.0 ;
	
	gettimeofday(&start_time, NULL);
	total_vertices=n_vertices;
	n_bcc = edgelists.size();
	
	total_bcc= n_bcc;
	newbcccnt= n_bcc;
	
	{
		
		for (int i = 0; i < n_bcc; i++ )
			bcc_sizes.pb(mp(edgelists[i].size(), bicc.bccmapping[i]));	//CH01
		sort(bcc_sizes.rbegin(), bcc_sizes.rend());

		vvii new_id_list = vvii(n_vertices);
		cout<<" n_bcc= "<<n_bcc<<endl;
		//exit(1);
		for (int i=0; i<n_bcc; i++) {
			int cnt = 0;
			map<int,int> new_id_map;
			map<int,pair<int,int> > new_old;
			map<int,int> rev_new_old;	//reverse new_old
			//vii new_old;
	
			for(int j = 0 ; j < edgelists[i].size(); j++ ){
				int u,v,w;
				u = edgelists[i][j].first.first;
				v = edgelists[i][j].first.second;
				w = edgelists[i][j].second;
	
				if(bicc.artmap.find(u)!=bicc.artmap.end())
				{
					if(bicc.bccmap.find(i)!=bicc.bccmap.end())					
					{
						(bicc.bccmap.find(i)->second).push_back(bicc.artmap.find(u)->second);					
					}
					if(bicc.bccmap.find(i)==bicc.bccmap.end())
					{
						vi temp_vector;
						temp_vector.push_back(bicc.artmap[u]);
						bicc.bccmap.insert(pair<int,vi>(i,temp_vector));
					}
					if(bicc.bccmap.find(bicc.artmap[u])!=bicc.bccmap.end())					
					{
						(bicc.bccmap.find(bicc.artmap[u])->second).push_back(i);					
					}
					if(bicc.bccmap.find(bicc.artmap[u])==bicc.bccmap.end())
					{
						vi temp_vector;
						temp_vector.push_back(i);
						bicc.bccmap.insert(pair<int,vi>(bicc.artmap[u],temp_vector));
					}	
			
				}
				
				if(bicc.artmap.find(v)!=bicc.artmap.end())
				{
					if(bicc.bccmap.find(i)!=bicc.bccmap.end())					
					{
						(bicc.bccmap.find(i)->second).push_back(bicc.artmap.find(v)->second);					
					}
					if(bicc.bccmap.find(i)==bicc.bccmap.end())
					{
						vi temp_vector;
						temp_vector.push_back(bicc.artmap[v]);
						bicc.bccmap.insert(pair<int,vi>(i,temp_vector));
					}
					if(bicc.bccmap.find(bicc.artmap[v])!=bicc.bccmap.end())					
					{
						(bicc.bccmap.find(bicc.artmap[v])->second).push_back(i);					
					}
					if(bicc.bccmap.find(bicc.artmap[v])==bicc.bccmap.end())
					{
						vi temp_vector;
						temp_vector.push_back(i);
						bicc.bccmap.insert(pair<int,vi>(bicc.artmap[v],temp_vector));
					}	
			
				}

				if(new_id_map.find(u) == new_id_map.end()) {
					new_id_map[u] = cnt;
					new_old[cnt]=mp(u,-1);
					if(bicc.artmap.find(u)!=bicc.artmap.end())
						rev_new_old[u]=cnt;
					bicc.block_index[u].second=cnt;
					//new_old.pb(mp(cnt,u));
					new_id_list[u].pb(mp(i,cnt));
					cnt++;
				}
				if(new_id_map.find(v) == new_id_map.end()) {
					new_id_map[v] = cnt;
					new_old[cnt]=mp(v,-1);
					if(bicc.artmap.find(v)!=bicc.artmap.end())
						rev_new_old[v]=cnt;
					bicc.block_index[v].second=cnt;
					//new_old.pb(mp(cnt,v));
					new_id_list[v].pb(mp(i,cnt));
					cnt++;
				}
				edgelists[i][j].first.first = new_id_map[u];
				edgelists[i][j].first.second = new_id_map[v];
			}
			bicc.bccvertno[i]=new_old;
			bicc.rev_bccvertno[i]=rev_new_old;
			bicc.bccvrtcnt.insert(mp(i,cnt));
		//	apsp_graphs.pb(APSP_Graph(cnt, edgelists[i]));
		}

		bicc.fillbccparent(graph.n_vertices);	
		bicc.reachfunction(graph.n_vertices,bc_result);
	
		bicc.to_eargraph(graph.n_vertices,graph.n_edges,edgelists,bicc.bccvertno,bicc.artmap,n_bcc, bicc.bccvrtcnt, bicc.inv_art, bc_result, g_free_time,
		active_vert_count, free_vert_count);


	}
}