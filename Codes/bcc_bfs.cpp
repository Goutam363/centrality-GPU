/*
bcc_bfs.cpp:582:17: error: expected ‘declare’ before ‘parallel’
 #pragma omp end parallel
                 ^
bcc_bfs.cpp:715:17: error: expected ‘declare’ before ‘parallel’
 #pragma omp end parallel
                 ^
bcc_bfs.cpp:854:17: error: expected ‘declare’ before ‘parallel’
 #pragma omp end parallel
*/

#include <iostream>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <list>
#include <ctime>
#include <map>
#include <omp.h>
#include <cstdlib>
#include "ear_graph.h"

#define PARALLEL_CUTOFF 10000
//#include <unordered_map>
#define MAX_UINT 999999999
#define MAX_LEVELS 1024

unsigned long long offset=4294967295;
unsigned long long bigoffset=18446744069414584320ULL;

using namespace std;

struct earVertex
{
	int leftDist;
	int postProcnode;
	int rightDist;
	unsigned int earWeight;

//	struct earVertex *nxtaddr;
};
struct adjList
{
	int count2D;
	unsigned long long endVertex;
	struct earVertex *nxtearVertex;
//	struct adjList *nxtadjList; 
};


void findnontree_edges(bool *&non_tree_e, int *&parents, int *&levels, int num_verts, int number_edges, int *&C, unsigned int *&R)
{

	int u,v,end,pos;
	for(int i=0;i<num_verts;i++)
	{

		u = parents[i];
		v = i;

	//	std::cout<<"u "<<u<<" v "<<v<<std::endl;

		if(u!=num_verts-1)
			end = R[u+1];
		else
			end = number_edges;
		for(int j=R[u]; j < end; j++)	
		{

			if(C[j]==v)
			{
				non_tree_e[j] = 1;

				if(v!=num_verts-1)
					pos = R[v+1];
				else
					pos = number_edges;

				for(int k=R[v]; k < pos ;k++)
				{
					if(C[k]==u)
					{
						non_tree_e[k] = 1;
					}	
				}	
			}	
		}	

	}	

}
void bcc_bfs_do_bfs(graph_ear_decompos& g, int root, int* parents, int* levels, int**& level_queues, int* level_counts, int& num_levels, 
vector<ear_premble*>&ear_premble_vec, vector<bool>& mark_eargraph_nodes, vector<bool>& mark_notfree_nodes);


clock_t t1,t2;
int earCount=0;; //Stores number of ears having 2 degree vertices 
unsigned int leftrightCount=0;
int preSum(list<unsigned long long>leftList,list<unsigned long long>rightList,unsigned int wtSum,unsigned int lWt,
unsigned int rWt,unsigned int nteWT,fstream &postProcout,int node1,int node2, vector<ear_premble*>&ear_premble_vec, 
vector<bool>& mark_eargraph_nodes ,vector<bool>& mark_notfree_nodes)
{
	
	unsigned int Weight,wt;
	wt=nteWT;
	int node;
	
	int index = 0;
	if(leftList.empty()==true && rightList.empty()==true)
	{

		return 0;
	}
	++earCount;
	leftrightCount+=(leftList.size()+rightList.size());
	postProcout <<-1<<" "<<node1<<" "<<node2<<" "<<(leftList.size()+rightList.size())<<" "<<wtSum<<endl;
	//std::cout<<-1<<" "<<node1<<" "<<node2<<" "<<(leftList.size()+rightList.size())<<" "<<wtSum<<std::endl;

	mark_eargraph_nodes[node1] = true;
	mark_eargraph_nodes[node2] = true;
	mark_notfree_nodes[node1] = true;
	mark_notfree_nodes[node2] = true;

	int *ld = new int[leftList.size() + rightList.size()];
	int *rd = new int[leftList.size() + rightList.size()];
	int *nodes = new int[leftList.size() + rightList.size()];
	//ear_p->lnode = node1;
	//ear_p->rnode = node2;
	//ear_p->cnt = leftList.size() + rightList.size();
	//ear_p->td = wtSum;
	ear_premble *ear_p = new ear_premble(node1, node2, leftList.size() + rightList.size(), wtSum);


	if(leftList.empty()!=true && rightList.empty()==true)
	{
		if(nteWT==0)
			wt=leftList.front()&offset;
		while(!leftList.empty())
		{
			Weight=leftList.front()&offset;
			node=(leftList.front() >> 32);
			//if(node==5)
//			postProcout<<" wt "<<wt<<" wtSum "<<wtSum<<" Weight "<<Weight<<" node "<<node<<endl;
			postProcout<<(wtSum - Weight)<<" "<<node<<" "<<wt<<endl;
		//	cout<<(wtSum - Weight)<<" "<<node<<" "<<wt<<endl;
			
			mark_notfree_nodes[node] = true;
			ear_p->ld[index] = (wtSum - Weight);
			ear_p->rd[index] = wt;
			ear_p->nodes[index] = node;
			index++;
			
			wt+=Weight;
			wtSum-=Weight;
			leftList.pop_front();
		} 
	}
	else if(leftList.empty()==true && rightList.empty()!=true)
	{	
		if(nteWT==0)
			wt=rightList.front()&offset;
		while(!rightList.empty())
		{
			Weight=rightList.front()&offset;
			node=(rightList.front() >> 32);
			
		//	postProcout<<" wt "<<wt<<" wtSum "<<wtSum<<" Weight "<<Weight<<" node "<<node<<endl; 
			postProcout<<wt<<" "<<node<<" "<<(wtSum - Weight)<<endl;
		//	cout<<wt<<" "<<node<<" "<<(wtSum - Weight)<<endl;
			
			mark_notfree_nodes[node] = true;
			ear_p->ld[index] = wt;
			ear_p->rd[index] = (wtSum - Weight);
			ear_p->nodes[index] = node;
			index++;
			
			wt+=Weight;
			wtSum-=Weight;
			rightList.pop_front();
		}

	}	
	else
	{
		wt=0;
		while(!leftList.empty())
		{
			Weight=leftList.front()&offset;
			node=(leftList.front() >> 32);
			

			postProcout<<(wtSum -(nteWT+rWt+wt))<<" "<<node<<" "<<(wt+nteWT+rWt)<<endl;
			//cout<<(wtSum -(nteWT+rWt+wt))<<" "<<node<<" "<<(wt+nteWT+rWt)<<endl;
			
			mark_notfree_nodes[node] = true;
			ear_p->ld[index] = (wtSum -(nteWT+rWt+wt));
			ear_p->rd[index] = (wt+nteWT+rWt);
			ear_p->nodes[index] = node;
			index++;

			wt+=Weight;
			//wtSum-=Weight;
			leftList.pop_front();
		}
		wt=0;
		while(!rightList.empty())
		{
			Weight=rightList.front()&offset;
			node=(rightList.front() >> 32);
		//	if(node==5)
		//		postProcout<<"wt "<<wt<<"nteWT "<<nteWT<<"lWt "<<lWt<<endl;
					
			postProcout<<(wt+nteWT+lWt)<<" "<<node<<" "<<(wtSum -(nteWT+lWt+wt))<<endl;
			//cout<<(wt+nteWT+lWt)<<" "<<node<<" "<<(wtSum -(nteWT+lWt+wt))<<endl;
		
			mark_notfree_nodes[node] = true;
			ear_p->ld[index] = (wt+nteWT+lWt);
			ear_p->rd[index] = (wtSum -(nteWT+lWt+wt));
			ear_p->nodes[index] = node;
			index++;

			wt+=Weight;
			//wtSum-=Weight;
			rightList.pop_front();
		}

	}

	//ear_p->ld = ld;
	//ear_p->rd = rd;
	//ear_p->nodes = nodes;
	ear_premble_vec.push_back(ear_p);
	//cout<<" Test "<<ear_premble_vec[0]->lnode<<" "<<ear_premble_vec[0]->rnode<<endl;
	return 0;	
}

string NumToStr(int num)
{
	stringstream ss;
	ss << num;
	string s=ss.str();
	return s;
}


void bcc_bfs(graph_ear_decompos& g, vector<ear_premble*>&ear_premble_vec, vector<bool>& mark_eargraph_nodes, vector<bool>& mark_notfree_nodes)
{

	int num_verts = g.num_vertices();
	int* parents = (int*)malloc(num_verts*sizeof(int));
	int* levels = (int*)malloc(num_verts*sizeof(int));
	int* high = (int*)malloc(num_verts*sizeof(int));
	int* low = (int*)malloc(num_verts*sizeof(int));
	bool* art_point = (bool*)malloc(num_verts*sizeof(bool));

	int num_levels = 0;
	int* level_counts = (int*)malloc(MAX_LEVELS*sizeof(int));
	int** level_queues = (int**)malloc(MAX_LEVELS*sizeof(int*));

	#pragma omp parallel for
	for (int i = 0; i < num_verts; ++i)
	{
		levels[i] = -1;
		parents[i] = -1;
		high[i] = -1;
		low[i] = -1;
		art_point[i] = false;
	}
	
	int root = g.get_max_degree_vert();

	levels[root] = 0;
	parents[root] = root;
	high[root] = root;

	bcc_bfs_do_bfs(g, root,	parents, levels, level_queues, level_counts, num_levels, ear_premble_vec, mark_eargraph_nodes, mark_notfree_nodes);

/*	printf("Level array \n");
	for(int x=0; x<num_verts; x++)
	{
		printf(" %d Level is %d\n",x,levels[x]);
	}

	printf("Parent array \n");
	for(int x=0; x<num_verts; x++)
	{
		printf(" %d parent is %d\n",x,parents[x]);
	}*/

	free(parents);
	free(levels);
	free(high);
	free(low);
	free(art_point);
}



void bcc_bfs_do_bfs(graph_ear_decompos& g, int root, int* parents, int* levels, int**& level_queues, int* level_counts, int& num_levels, 
vector<ear_premble*>&ear_premble_vec, vector<bool>& mark_eargraph_nodes, vector<bool>& mark_notfree_nodes)
{
	//cout<<"Time "<<ctime(&time(0))<<endl;
	
	int num_verts = g.num_vertices();
	double avg_out_degree = g.get_avg_out_degree();
	bool already_switched = false;

	int num_threads = omp_get_max_threads();
	cout<<"Num max threads "<<num_threads<<endl;
	int* to_visit = new int[num_verts];
	int back = 0;	
	int num_descs = 0;
	int* to_visit_threads = new int[num_threads*num_verts];
	int* thread_backs = new int[num_threads];

	bool* visited = new bool[num_verts];
	for (int i = 0; i < num_verts; ++i)
		visited[i] = false;
	
	double elt, elt2;
	double alpha = 15.0;
	double beta = 25.0;

	int v = root;
	parents[v] = v;
	levels[v] = 0;
	level_queues[0] = (int*)malloc(1*sizeof(int));
	level_queues[0][0] = v;
	level_counts[0] = 1;

	to_visit[0] = v;
	back = 1;
	bool *non_tree_e=g.get_out_non_tree_edge();				//A001
	unsigned int *out_degree_l=g.get_out_degree_list(); 	//A001
	unsigned int *master=g.get_master(); 					//A001
	unsigned int *serial=g.get_serial(); 					//A001
	unsigned int *lca_array=g.get_lca_array();
	int *out_arr=g.get_out_vertices();						//A001
	int *weights_arr=g.get_out_weights();					//A001
	unsigned int pos,pos1; 										//A001
	int number_edges=g.num_edges();							//A001
	unsigned int limits,limits1,i,j,k,j1;									//A001
//	struct find_master *s_master=(struct find_master *)malloc(sizeof(struct find_master)*number_edges);
		
	int local_num_descs;
	bool use_hybrid = false;

t1=clock();	
double omptime=omp_get_wtime();
//omp_set_num_threads(1);
#pragma omp parallel 
{
	int tid = omp_get_thread_num();		
	int* to_visit_thread = &to_visit_threads[tid*num_verts];
	int thread_back = 0;
	int level = 1;
	int prev_level;
	int vert, out_degree, out;
	int* outs;
	while (back)
	{
#if DEBUG
		elt = timer();
#endif
		local_num_descs = 0;
		thread_back = 0;

		if (!use_hybrid)
		{
if (num_verts - num_descs > PARALLEL_CUTOFF)
{
	#pragma omp for schedule(static) reduction(+:local_num_descs)
			for (int i = 0; i < back; ++i)
			{
				vert = to_visit[i];
				out_degree = g.out_degree(vert);
				outs = g.out_vertices(vert);
				for (int j = 0; j < out_degree; ++j)
				{			
					out = outs[j];	
					if (levels[out] < 0)
					{
						levels[out] = level;
						parents[out] = vert;
					//	printf(" parent=%d child=%d\n",vert,out );	//A001
					//	printf("vertex=%d \n",vert);
						to_visit_thread[thread_back++] = out;
						++local_num_descs;
						
						pos=out_degree_l[vert];						//A001
				
						if(vert!=(num_verts-1))						//A001
							limits=out_degree_l[vert+1];			//A001
						else										//A001
							limits=number_edges;					//A001

						for(unsigned int k1=pos;k1<limits;k1++)		//A001
						{				
							if(out_arr[k1]==out)					//A001
								non_tree_e[k1]=1;					//A001 1 indicate tree edge
						}

						//For B-->A	
						pos=out_degree_l[out];						//A001
						if(out!=(num_verts-1))						//A001
							limits=out_degree_l[out+1];				//A001
						else										//A001
							limits=number_edges;					//A001	

						for(unsigned int k1=pos;k1<limits;k1++)		//A001
						{				
							if(out_arr[k1]==vert)					//A001
								non_tree_e[k1]=1;					//A001 1 indicate tree edge
						}	
					}
				}
			}
}
else
{
	#pragma omp single
	{
			for (int i = 0; i < back; ++i)
			{
				vert = to_visit[i];
				out_degree = g.out_degree(vert);
				outs = g.out_vertices(vert);
				for (int j = 0; j < out_degree; ++j)
				{			
					out = outs[j];	
					if (levels[out] < 0)
					{
						levels[out] = level;
						parents[out] = vert;
						to_visit_thread[thread_back++] = out;
						++local_num_descs;
					//	printf(" parent=%d child=%d\n",vert,out );	//A001
						//For A-->B
						//To mark parent to child edge
						pos=out_degree_l[vert];						//A001
						
						if(vert!=(num_verts-1))						//A001
							limits=out_degree_l[vert+1];			//A001	//To find out range of destination vertices
						else										//A001
							limits=number_edges;					//A001	//For thr last vertex

						for(unsigned int k1=pos;k1<limits;k1++)		//A001
						{				
							if(out_arr[k1]==out)					//A001 
								non_tree_e[k1]=1;					//A001 1 indicate tree edge 
						}
						//For B-->A		For undirected graph A to B and B to A both edges are present so here we are marking edge B to A
						//To mark child to parent edge
						pos=out_degree_l[out];						//A001
						if(out!=(num_verts-1))						//A001
							limits=out_degree_l[out+1];				//A001
						else										//A001
							limits=number_edges;					//A001	


						for(unsigned int k1=pos;k1<limits;k1++)		//A001
						{				
							if(out_arr[k1]==vert)					//A001 
								non_tree_e[k1]=1;					//A001 1 indicate tree edge
						}

					}
				}
			}
	}
}
		}
		else
		{
			prev_level = level - 1;

	#pragma omp for schedule(static) reduction(+:local_num_descs)
			for (int i = 0; i < num_verts; ++i)
			{
				vert = i;
				//printf("levels[vert]=%d vert=%d\n",levels[vert],vert);
				if (levels[vert] < 0)
				{//	printf("vert=%d\n",vert);				//A001
					out_degree = g.out_degree(vert);
					outs = g.out_vertices(vert);
					for (int j = 0; j < out_degree; ++j)
					{
						out = outs[j];
						if (levels[out] == prev_level)
						{
							levels[vert] = level;
							parents[vert] = out;
							to_visit_thread[thread_back++] = vert;
							++local_num_descs;
					//		printf(" parent=%d child=%d\n",vert,out );	//A001
							pos=out_degree_l[vert];						//A001
				
							if(vert!=(num_verts-1))						//A001
								limits=out_degree_l[vert+1];			//A001
							else										//A001
								limits=number_edges;					//A001
						
							for(unsigned int k1=pos;k1<limits;k1++)		//A001
							{				
								if(out_arr[k1]==out)					//A001
									non_tree_e[k1]=1;					//A001 1 indicate tree edge
							}
						
							//For B-->A	
							pos=out_degree_l[out];						//A001
				
							if(out!=(num_verts-1))						//A001
								limits=out_degree_l[out+1];				//A001
							else										//A001
								limits=number_edges;					//A001	

							for(unsigned int k1=pos;k1<limits;k1++)		//A001
							{				
								if(out_arr[k1]==vert)					//A001
									non_tree_e[k1]=1;					//A001 1 indicate tree edge
							}	
							break;
						}
					}
				}
			}
		}
		thread_backs[tid] = thread_back;

#pragma omp barrier


if (tid == 0)
{		
		num_descs += local_num_descs;

#if DEBUG
		printf("num_descs: %d   local: %d\n", num_descs, local_num_descs);
#endif

		if (!use_hybrid)
		{	
			double edges_frontier = (double)local_num_descs * avg_out_degree;
			double edges_remainder = (double)(num_verts - num_descs) * avg_out_degree;
			if ((edges_remainder / alpha) < edges_frontier && edges_remainder > 0 && !already_switched)
			{
#if DEBUG	 
				printf("\n=======switching to hybrid\n\n");
#endif
				use_hybrid = true;
			}
#if DEBUG
			printf("edge_front: %f  edge_rem: %f\n", edges_frontier, edges_remainder);
#endif
		}
		else
		{
			if ( ((double)num_verts / beta) > local_num_descs  && !already_switched)
			{
#if DEBUG
				printf("\n=======switching back\n\n");
#endif
				use_hybrid = false;
				already_switched = true;
			}
		}
		
#if DEBUG
		elt2 = timer();
#endif
		back = 0;
		for (int i = 0; i < num_threads; ++i)
		{
			if (thread_backs[i])
			{
				int offset = i*num_verts;
				for (int j = 0; j < thread_backs[i]; ++j)
				{	to_visit[back++] = to_visit_threads[offset + j];
					//printf("to visit %d\n",to_visit_threads[offset + j] );
				}	
			}
		}

		level_queues[level] = (int*)malloc(back*sizeof(int));
		copy(to_visit, to_visit+back, level_queues[level]);
		level_counts[level] = back;
		++num_levels;
#if DEBUG
		elt2 = timer() - elt2;
		printf("create array: %9.6lf\n", elt2);
#endif
	
#if DEBUG
		elt = timer() - elt;
		printf("level: %d  num: %d  time: %9.6lf\n", level, back, elt);	
#endif
}
	++level;

#pragma omp barrier

	}
} // end parallel
//#pragma omp end parallel

omptime=omp_get_wtime()-omptime;
t2=clock();
//printf("Time required for before LCA step is %f %lf\n", ((float(t2)-float(t1))/CLOCKS_PER_SEC),omptime);

findnontree_edges(non_tree_e,parents,levels,num_verts,number_edges, out_arr, out_degree_l);
//cout<<"Time NTE "<<((float(t2)-float(t1))/CLOCKS_PER_SEC)<<endl;
//=======
	for(i=0;i<num_verts;i++)
	{
		pos=out_degree_l[i];
		if(i!=(num_verts-1))
			limits=out_degree_l[i+1];
		else
			limits=number_edges;

		//printf("i=%d threadId=%d\n",i,omp_get_thread_num());
		for(j=pos;j<limits;j++)
		{
			if(non_tree_e[j]==0)					//To check If it is non-tree edge 
			{

				pos1=out_degree_l[out_arr[j]];
				if(j!=(num_verts-1))
					limits1=out_degree_l[out_arr[j]+1];
				else
					limits1=number_edges;
				
				for(j1=pos1;j1<limits1;j1++)
				{	
					if(out_arr[j1]==i)
						non_tree_e[j1]=1;			
				}

				if(parents[i]==out_arr[j] || parents[out_arr[j]]==i)
				{
					printf(" u %d parents[%d] = %d v %d parents[%d] = %d level diff %d\n",i,i,parents[i], out_arr[j], out_arr[j], parents[out_arr[j]], abs( levels[i] - levels[out_arr[j]] ) );
					exit(0);
				}

			}	
		}
	}
//=======
#if CDEBUG	
for (i = 0; i < num_verts; i++)
         printf(" o_l=%u \n",out_degree_l[i]);
#endif
// exit(0);
//==============FIND LCA()===========// A001
//***********struct get_lca *lca=(struct get_lca *)malloc(sizeof(struct get_lca)*number_edges);
set<unsigned long long>ntEdges; //Set to store concate result of edge(u,v) i.e. uv 
set<unsigned long long>::iterator it;  //Try to use unorder set it may improve performance
unsigned long long strUV,strVU;
unsigned int l2=0,dh,k1,l3;//i,j,j1;
static int l1=-1;
int u_node,v_node,h1,h2,base1,base2; 
static int threadId=-1;
int dnc_v_u_edge;  //A001    //Do not consider v to u edge bcz we have already considered u to v edge
	
	//omp_set_dynamic(1);     // Explicitly disable dynamic teams
	t1=clock();
	omptime=omp_get_wtime();	
	omp_set_num_threads(1);
	//printf(" %d %d\n",l1,omp_get_thread_num()); 
	#pragma omp parallel for private(j,pos,limits,u_node,v_node,h1,h2,base1,base2,strUV,strVU,it,dnc_v_u_edge,dh,threadId,k1) schedule(static)//lastprivate(l1)		
	for(i=0;i<num_verts;i++)
	{
		pos=out_degree_l[i];
		if(i!=(num_verts-1))
			limits=out_degree_l[i+1];
		else
			limits=number_edges;

		//printf("i=%d threadId=%d\n",i,omp_get_thread_num());
		for(j=pos;j<limits;j++)
		{
		//	printf("edge %d to %d = %d\n",i,out_arr[j],non_tree_e[j]);
			if(non_tree_e[j]==0)					//To check If it is non-tree edge 
			{
				 u_node=i;
				 v_node=out_arr[j];	
				 h1=levels[u_node];
				 h2=levels[v_node];
				 base1=u_node;
				 base2=v_node;

				if(parents[u_node]==v_node || parents[v_node]==u_node)	
				{
					printf("=====u_node %d ===== v_node %d ==== lev %d\n",u_node,v_node,abs(h1-h2));
				}	

				if(h1 > h2)
			 	{
			
					 dh=h1-h2;
					for(k1=0; k1<dh; k1++)
					{
						//base2 = v;
						u_node = parents[u_node];
					}	
				}
				if(h2>h1)
				{
					 dh=h2-h1;
					for(k1=0; k1<dh; k1++)
					{
						//base2 = v;
						v_node = parents[v_node];
					}	
				}
				
				while(1)
				{
					if(u_node==v_node)
					{
						
						//if(dnc_v_u_edge==0)		//Do not consider v to u edge if u to v edge already inserted
							lca_array[j]=u_node;
							//printf("u=%u v=%u LCA=%u\n",i,out_arr[j],u_node);		
							break;
		
					}else
					{
						v_node=parents[v_node];
						u_node=parents[u_node];
					}
		
				}	
			}
		}		
	}
//#pragma omp end parallel
l1=l1+1;		//Because in above code l1 is increamented befor adding data to lca data structure
//printf("hi4==%d\n",l1);
t2=clock();
omptime=omp_get_wtime()-omptime;
//cout<<"Time LCA "<<((float(t2)-float(t1))/CLOCKS_PER_SEC)<<endl;
//printf("Time required for Find LCA step is %f \n", ((float(t2)-float(t1))/CLOCKS_PER_SEC));
//printf("Time required for Find LCA step is using Openmp Timer %lf\n",omptime);

ntEdges.clear();
//==============FIND LCA() ENDS===========// A001 	A001

/*//=====For test=====
int k2=l1;
for(i=0;i<k2;i++)
{
	printf("u=%d v=%d lca=%d ser=%u w=%d\n",lca[i].u1,lca[i].v1,lca[i].lca_u1_v1,lca[i].serial_number,lca[i].e_weight);
}
//=====End for test=====*/

//==============Assign master to every tree-edge=========

//printf("hi7 %d\n",s_master[0].master_arr[0]);
int ntv1,ntv2;
omptime=omp_get_wtime();
int u1_node,v1_node,lca_node,level_lca;
//omp_set_num_threads(1);
//#pragma omp parallel for //private(j,j1)
for(i=0;i<num_verts;i++)
{
	
	pos=out_degree_l[i];
	if(i!=(num_verts-1))
		limits1=out_degree_l[i+1];
	else
		limits1=number_edges;
	for(k=pos;k<limits1;k++)
	{
		if(non_tree_e[k]==0)					//To check If it is non-tree edge 
		{
			u1_node=i;
			v1_node=out_arr[k];
			lca_node=lca_array[k];
			level_lca=levels[lca_node];
	
		//	printf("u1_node %d v1_node %d lca %d parent[u1_node] %d parent[v1_node] %d\n",u1_node,v1_node,lca_node, parents[u1_node],parents[v1_node]);	
			while(1)
			{
		//		printf("Test1\n");
				int parent_u1=parents[u1_node];
				pos=out_degree_l[parent_u1];
				if(parent_u1!=(num_verts-1))						//A001
					limits=out_degree_l[parent_u1+1];				//A001
				else												//A001
					limits=number_edges;							//A001
								
				for(j=pos;j<limits;j++)
				{
					if((out_arr[j]==u1_node) && (!(non_tree_e[j]==0)) )
					{
						if(master[j]==999999999)
						{
							master[j]=k;			//Number of non-tree edge
						//	if(u1_node==237 && parent_u1==15457)
						//		printf("master %d %d edge %d %d %d %u lca %d hi5\n",i,out_arr[k],parent_u1,u1_node,k,serial[k],lca_node);
						//	if(u1_node==15457 && parent_u1==228)
						//		printf("master %d %d edge %d %d %d %u lca %d hi5\n",i,out_arr[k],parent_u1,u1_node,k,serial[k],lca_node);
		//					printf("Test11\n");
						}
						else
						{
						//	printf("master=%u %u i=%u j=%u %u\n",lca_array[master[j]],master[j],i,k,out_arr[j]);	
						//	if(u1_node==237 && parent_u1==15457)
						//				printf("master %d %d edge %d %d %d %u lca %d hi4\n",i,out_arr[k],parent_u1,u1_node,k,serial[k],lca_node);
						//	if(u1_node==15457 && parent_u1==228)
						//				printf("master %d %d edge %d %d %d %u lca %d hi4\n",i,out_arr[k],parent_u1,u1_node,k,serial[k],lca_node);		

							if(levels[lca_array[master[j]]] > level_lca)
							{
						//		printf("Test12\n");
								master[j]=k;
						//		if(u1_node==237 && parent_u1==15457)
						//			printf("master %d %d edge %d %d %d lca %d hi2 \n",i,out_arr[k],parent_u1,u1_node,k,lca_node);
							}
							else
								if((levels[lca_array[master[j]]] == level_lca) && ( serial[master[j]] > serial[k]))
								{
									master[j]=k;	
						//			if(u1_node==237 && parent_u1==15457)
						//				printf("master %d %d edge %d %d %d %u lca %d hi3\n",i,out_arr[i],parent_u1,u1_node,k,serial[k],lca_node);
								}
		//					
						}

					//printf("\n u= %d v=%d m=%d threadId=%d s_master=%d",s_master[u1_node].unode,u1_node,i,omp_get_thread_num(),s_master[u1_node].master_arr[i]);
					}
				}
				if(parent_u1==lca_node)
					break;
				u1_node=parent_u1;
			}

			while(1)
			{
		//			printf("Test2\n");
				int parent_v1=parents[v1_node];
				pos=out_degree_l[parent_v1];
				if(parent_v1!=(num_verts-1))						//A001
					limits=out_degree_l[parent_v1+1];				//A001
				else												//A001
					limits=number_edges;							//A001
								
				for(j=pos;j<limits;j++)
				{
					if((out_arr[j]==v1_node) && (!(non_tree_e[j]==0)) )
					{
						if(master[j]==999999999)
						{
							master[j]=k;			//Number of non-tree edge
						//	if(v1_node==237 && parent_v1==15457)
						//		printf("master %d %d edge %d %d %d %u lca %d hi55\n",i,out_arr[k],parent_v1,v1_node,k,serial[k],lca_node);
						//	if(v1_node==15457 && parent_v1==228)
						//		printf("master %d %d edge %d %d %d %u lca %d hi55\n",i,out_arr[k],parent_v1,v1_node,k,serial[k],lca_node);
						}
						else
						{
						//	printf("master22=%u\n",master[j]);	
						//	if(v1_node==237 && parent_v1==15457)
						//				printf("master %d %d edge %d %d %d %u lca %d hi44\n",i,out_arr[k],parent_v1,v1_node,k,serial[k],lca_node);
						//	if(v1_node==15457 && parent_v1==228)
						//				printf("master %d %d edge %d %d %d %u lca %d hi44\n",i,out_arr[k],parent_v1,v1_node,k,serial[k],lca_node);		

							if(levels[lca_array[master[j]]] > level_lca)
							{
								master[j]=k;
						//		if(v1_node==237 && parent_v1==15457)	
						//			printf("master %d %d edge %d %d %d lca %d hi22\n",i,out_arr[k],parent_v1,v1_node,k,lca_node);
							}
							else
								if((levels[lca_array[master[j]]] == level_lca) && ( serial[master[j]] > serial[k]))
								{	
									master[j]=k;	
						//			if(v1_node==237 && parent_v1==15457)
						//				printf("master %d %d edge %d %d %d %u lca %d hi33\n",i,out_arr[k],parent_v1,v1_node,k,serial[k],lca_node);
								}
						}

					//	printf("\nu= %d v=%d m=%d threadId=%d s_master=%d",s_master[v1_node].unode,v1_node,i,omp_get_thread_num(),s_master[v1_node].master_arr[i]);	
					}
				}
				if(parent_v1==lca_node)
					break;
				v1_node=parent_v1;
			}
		}
	}	
}


//printf("hi8\n");
//#pragma omp end parallel

//==============Master ENDS===========// A001 	A001
omptime=omp_get_wtime()-omptime;
//printf("Time required for Assigned Master step using omptimer is %lf \n", omptime);

//==============Print edges in ear========= 	A001

//int u1_node,v1_node,lca_node,level_lca;
//hash<string,unsigned int>vertDistMap;						//To store distance of vertex from end points of ear
fstream postProcout,fout,foutear;
postProcout.open("PostProcessOp.txt",ios::out);
fout.open("Output.txt",ios::out);
list<unsigned long long>leftvertWtList;
list<unsigned long long>leftvertWtList1;
list<unsigned long long>rightvertWtList;
//vector<ear_premble> ear_premble_vec;
string strVertex;	
unsigned long long postprocVrtWt;
unsigned int weight_sum=0,udegreeWeightSum;
unsigned int wtLeft=0,wtRight=0,nteWT=0;
unsigned int postprocWt=0;
//int store_u1_ear,store_v1_ear;//break_flag1,break_flag2;
int save_u1_node,saveFirst_u1_node,parent_u1,parent_v1,tid;
bool print_flag;
int temp_parent_v1=MAX_UINT;
//cout<<"hi111"<<endl;
//omp_set_num_threads(1);
//#pragma omp parallel for private(j,udegreeWeightSum,save_u1_node,saveFirst_u1_node,print_flag,u1_node,v1_node,lca_node,level_lca,weight_sum,parent_u1,parent_v1,pos,limits,strVertex,postprocVrtWt,wtLeft,wtRight,nteWT,temp_parent_v1) schedule(static)

for(i=0;i<num_verts;i++)
{
	
	pos=out_degree_l[i];
	if(i!=(num_verts-1))
		limits1=out_degree_l[i+1];
	else
		limits1=number_edges;
	for(k=pos;k<limits1;k++)
	{
	//	if(i==0)
	//		printf("k=%u %d\n",k,limits1);
		if(non_tree_e[k]==0)					//To check If it is non-tree edge 
		{
			u1_node=i;
			v1_node=out_arr[k];
			ntv1 = u1_node;
			ntv2 = v1_node;
			lca_node=lca_array[k];
			level_lca=levels[lca_node];

			udegreeWeightSum=0;
			save_u1_node=MAX_UINT;
			saveFirst_u1_node=MAX_UINT;
			print_flag=false;
						
	//		if(i==0 && out_arr[k]==4)
	//			printf("LCA=%u\n",lca_node);
			weight_sum=0;
	//To process left side edges of non tree edge
			while(1)
			{
				parent_u1=parents[u1_node];
				pos=out_degree_l[parent_u1];
				if(parent_u1!=(num_verts-1))						//A001
					limits=out_degree_l[parent_u1+1];				//A001
				else												//A001
					limits=number_edges;							//A001
								
				for(j=pos;j<limits;j++)
				{
					if((out_arr[j]==u1_node) && (!(non_tree_e[j]==0)) )
					{
						if(ntv1==237 && ntv2==10303)
						{
							fout<<" left master "<<master[j]<<" k "<<k<<" parent_u1 "<<parent_u1<<" child "<<u1_node<<endl;
						}			

						if(master[j]==k)
						{
							//degree of u is greater than 2   	u---------------v non tree edge
							//This if and its else part can be combine 
							/*if(ntv1==237 && ntv2==10303)
							{
								fout<<" parent_u1 "<<parent_u1<<" child "<<u1_node<<" g.out_degree(i) "<<g.out_degree(i)<<endl;
							}*/	
							if(g.out_degree(i)>2)
							{
								weight_sum+=weights_arr[j];
								
								if(g.out_degree(u1_node)>2)
								{
									//if(save_u1_node!=MAX_UINT)
									//e(u,v) is non tree edge, saveFirst_u1_node stores first degree > 2 vertex from (u,l) path l is lca(u,v)
									if(saveFirst_u1_node==MAX_UINT)
										saveFirst_u1_node=u1_node;
											
									save_u1_node=u1_node;			//Initialize to MAX_UINT
									print_flag=true;
								}
								else
								{
									//For post processing part, storing vertex of degree 2 and its weight in map
									 postprocVrtWt= offset & u1_node;
									 postprocVrtWt = (postprocVrtWt << 32);
									 postprocVrtWt |= (offset & weights_arr[j]);
									 
									 leftvertWtList.push_back(postprocVrtWt);
								} 

								if(g.out_degree(parent_u1) > 2 && save_u1_node!=MAX_UINT)	
								{																		//	3	  2	   3
									//#pragma omp critical													y-----x----u---------v	
										//cout<<save_u1_node<<" "<<parent_u1<<" "<<weight_sum<<endl;	
									//printf("%d %d %u\n",save_u1_node,parent_u1,weight_sum);
									fout<<save_u1_node<<" "<<parent_u1<<" "<<weight_sum<<endl;
								/*	if(weight_sum==1 && (leftvertWtList.empty()!=true || rightvertWtList.empty()!=true))
									{
										fout<<" a ntv1 "<<ntv1<<" ntv2 "<<ntv2<<endl;	
									}
									if(ntv1==237 && ntv2==10303)
									{
										fout<<"u1 node "<<u1_node<<endl;
									}	*/
									preSum(leftvertWtList, rightvertWtList, weight_sum, 0, 0, nteWT, postProcout, save_u1_node, parent_u1, ear_premble_vec,
									mark_eargraph_nodes , mark_notfree_nodes);
									leftvertWtList.clear();
									rightvertWtList.clear();
									//reset weight_sum
									print_flag=false;
									weight_sum=0;
								}
							}
							//degree of u is less than 2   				u---------------v non tree edge
							else
							{
								weight_sum+=weights_arr[j];

								if(g.out_degree(u1_node)>2)
								{
									//if(save_u1_node!=MAX_UINT)
									//e(u,v) is non tree edge, saveFirst_u1_node stores first degree > 2 vertex from (u,l) path l is lca(u,v)
									if(saveFirst_u1_node==MAX_UINT)
										saveFirst_u1_node=u1_node;

									save_u1_node=u1_node;
									print_flag=true;
								}
								else
								{
									//For post processing part, storing vertex of degree 2 and its weight in map
									 postprocVrtWt = offset & u1_node;
									 postprocVrtWt = postprocVrtWt << 32;
									 postprocVrtWt |= (offset & weights_arr[j]);
									 
									 if(saveFirst_u1_node==MAX_UINT)
									 	leftvertWtList.push_back(postprocVrtWt);
									 else
									 	leftvertWtList1.push_back(postprocVrtWt);
									}
								
								if(save_u1_node==MAX_UINT)		//if(g.out_degree(lca[i].u1)>2 && save_u1_node==MAX_UINT)
								{
									udegreeWeightSum=weight_sum;	//Intialize to 0
								}

								if(g.out_degree(parent_u1) > 2 && save_u1_node!=MAX_UINT)	
								{
									//#pragma omp critical 
									//	cout<<save_u1_node<<" "<<parent_u1<<" "<<weight_sum<<endl;	
									
									//Instead of cout printf is used because in cout threads overlap each others data  
		//fout						printf("%d %d %u\n",save_u1_node,parent_u1,weight_sum);
									fout<<save_u1_node<<" "<<parent_u1<<" "<<weight_sum<<endl;						
									/*if(weight_sum==1 && (leftvertWtList1.empty()!=true || rightvertWtList.empty()!=true))
									{
										fout<<" b ntv1 "<<ntv1<<" ntv2 "<<ntv2<<endl;	
									}*/
									preSum(leftvertWtList1, rightvertWtList, weight_sum, 0, 0, nteWT, postProcout, save_u1_node, parent_u1, ear_premble_vec, 
									mark_eargraph_nodes , mark_notfree_nodes);
									leftvertWtList1.clear();
									rightvertWtList.clear();
									//reset weight_sum
									print_flag=false;
									weight_sum=0;
								}	
							}		
											
							
							//printf("%d %d\n",parent_u1,u1_node);			//Number of non-tree edge
						}
						else
						{

							//e(u,v) is non tree edge, saveFirst_u1_node stores first degree > 2 vertex from (u,l) path l is lca(u,v)
							if(saveFirst_u1_node==MAX_UINT)
								saveFirst_u1_node=u1_node;

							//printf("b1 %d %d m= %d %d master %d %d\n",store_u1_ear,store_v1_ear,lca[l1].u1,lca[l1].v1,master[j],i);
						}

					}
				}
				if(parent_u1==lca_node)
				{
					//printf("%d ",store_u1_ear);
					if(saveFirst_u1_node==MAX_UINT)
					{
						saveFirst_u1_node=parent_u1;	//To handle 2-4-7 triangle condition when NTE is (7,4) in testfile with V=9 E=22
					}

					break;
				}
				u1_node=parent_u1;
			}
		//	cout<<"Test2 "<<saveFirst_u1_node<<" "<<udegreeWeightSum<<" "<<lca[i].e_weight<<endl;
			//udegreeWeightSum+=lca[i].e_weight;
			wtLeft=udegreeWeightSum;
			weight_sum=weights_arr[k];//0;
			save_u1_node=MAX_UINT;
			//To process Right side edges of non tree edge
			while(1)
			{
				parent_v1=parents[v1_node];
				pos=out_degree_l[parent_v1];
				if(parent_v1!=(num_verts-1))						//A001
					limits=out_degree_l[parent_v1+1];				//A001
				else												//A001
					limits=number_edges;							//A001
								
				for(j=pos;j<limits;j++)
				{
					if((out_arr[j]==v1_node) && (!(non_tree_e[j]==0)) )
					{
						/*if(ntv1==237 && ntv2==10303)
						{
							fout<<" right master "<<master[j]<<" k "<<k<<" parent_v1 "<<parent_v1<<" child "<<v1_node<<endl;
						}*/	
						if(master[j]==k)
						{
							//degree of v is greater than 2
							/*if(ntv1==237 && ntv2==10303)
							{
								fout<<" parent_v1 "<<parent_v1<<" child "<<v1_node<<" g.out_degree(out_arr[k]) "<<g.out_degree(out_arr[k])<<endl;
							}*/
							if(g.out_degree(out_arr[k])>2)
							{
								if(saveFirst_u1_node==MAX_UINT)
									weight_sum+=weights_arr[j];
								
								if(g.out_degree(v1_node)>2)
								{
			//						if(saveFirst_u1_node==0)		//For testing purpose
			//							cout<<"Test1"<<saveFirst_u1_node<<endl; //For testing purpose
									//if(save_u1_node!=MAX_UINT)
									//if(lca[i].v1==u1_node)
															
									save_u1_node=v1_node;			//Initialize to MAX_UINT
									//This is to print first degree > 2 vertex from u to l and first degree > 2 vertex v to l
									if (saveFirst_u1_node!=MAX_UINT)
									{
										nteWT=weight_sum;
										save_u1_node=saveFirst_u1_node;
										weight_sum+=udegreeWeightSum;
										udegreeWeightSum=0;
										parent_v1=v1_node;

									}	
									
									saveFirst_u1_node=MAX_UINT;	
									
									print_flag=true;
								}
								else
								{
									//For post processing part, storing vertex of degree 2 and its weight in map
									 //nteWT=0;
									 postprocVrtWt = offset & v1_node;
									 postprocVrtWt = postprocVrtWt << 32;
									 //postprocVrtWt = bigoffset & weights_arr[j];
									 postprocVrtWt |= (offset & weights_arr[j]);
									 rightvertWtList.push_back(postprocVrtWt);

									//strVertex=NumToStr(v1_node);	
									//strVertex.append('R');		//Indicate Right side of NTE (non tree edge)
									//vertDistMap.insert(pair<string,int>(strVertex,weights_arr[j]));
								}

								if(g.out_degree(parent_v1) > 2 && save_u1_node!=MAX_UINT)	
								{
									//cout<<save_u1_node<<" "<<parent_v1<<" "<<weight_sum<<endl;	
		//fout						printf("%d %d %u\n",save_u1_node,parent_v1,weight_sum);
									fout<<save_u1_node<<" "<<parent_v1<<" "<<weight_sum<<endl;	
								/*	if(weight_sum==1 && (leftvertWtList.empty()!=true || rightvertWtList.empty()!=true))
									{
										fout<<" c ntv1 "<<ntv1<<" ntv2 "<<ntv2<<" u1 node "<<u1_node<<endl;	
									}*/
									//vertDistMap.insert(pair<string,int>("NTE",lca[i].e_weight));
									preSum(leftvertWtList, rightvertWtList, weight_sum, 0, 0, nteWT, postProcout, save_u1_node, parent_v1, ear_premble_vec, 
									mark_eargraph_nodes ,mark_notfree_nodes);
									rightvertWtList.clear();
									leftvertWtList.clear();
									nteWT=0;
									//reset weight_sum
									print_flag=false;
									weight_sum=0;
								}
							}
							//degree of v is less than 2
							else
							{
								if(saveFirst_u1_node==MAX_UINT)
									weight_sum+=weights_arr[j];
					
								if(g.out_degree(v1_node)>2)
								{
									//if(save_u1_node!=MAX_UINT)
									//	saveFirst_u1_node=v1_node;

									save_u1_node=v1_node;
									print_flag=true;
								}
								else
								{
									//For post processing part, storing vertex of degree 2 and its weight in map
									 
									 postprocVrtWt = offset & v1_node;
									 postprocVrtWt = postprocVrtWt << 32;
									 postprocVrtWt |= (offset & weights_arr[j]);
									 rightvertWtList.push_back(postprocVrtWt);

								}
								
								//This is to print first degree > 2 vertex from u to l and first degree > 2 vertex v to l
								if(saveFirst_u1_node!=MAX_UINT)
								{
									nteWT=weight_sum;
									save_u1_node=saveFirst_u1_node;
									weight_sum+=udegreeWeightSum;
									udegreeWeightSum=0;
									weight_sum+=weights_arr[j];
									saveFirst_u1_node=MAX_UINT;
									//temp_parent_v1=parent_v1;
									//parent_v1=v1_node;
								}
								//Check dis step
									saveFirst_u1_node=MAX_UINT;	


								if(g.out_degree(parent_v1) > 2 && save_u1_node!=MAX_UINT)	
								{
									//cout<<save_u1_node<<" "<<parent_v1<<" "<<weight_sum<<endl;	
		//fout						printf("%d %d %u\n",save_u1_node,parent_v1,weight_sum);
									fout<<save_u1_node<<" "<<parent_v1<<" "<<weight_sum<<endl;						
								/*	if(weight_sum==1 && (leftvertWtList.empty()!=true || rightvertWtList.empty()!=true))
									{
										fout<<" d ntv1 "<<ntv1<<" ntv2 "<<ntv2<<endl;	
									}*/
									//To avoid repeating of v node when its degree is < 2 i.e. when parent_v1=v1_node happen 
									wtRight=weight_sum-(wtLeft+weights_arr[k]);
									preSum(leftvertWtList, rightvertWtList, weight_sum, wtLeft, wtRight, nteWT, postProcout, save_u1_node, parent_v1, ear_premble_vec, 
									mark_eargraph_nodes ,mark_notfree_nodes);
									rightvertWtList.clear();
									leftvertWtList.clear();
									nteWT=0;
									print_flag=false;
									weight_sum=0;
								}	
							}		
							
						}
						else
						{
							if(saveFirst_u1_node==MAX_UINT)
								weight_sum+=weights_arr[j];
							
							if(saveFirst_u1_node!=MAX_UINT)
								{
									nteWT=weight_sum;
									if(udegreeWeightSum!=0)
										weight_sum+=udegreeWeightSum;
								//if(saveFirst_u1_node==10 && v1_node==12)
								//cout<<"weight= "<<weight_sum<<endl;	
										
		//fout						printf("%d %d %u\n",saveFirst_u1_node,v1_node,weight_sum);
									fout<<saveFirst_u1_node<<" "<<v1_node<<" "<<weight_sum<<endl;					
		//							preSum(rightvertWtList,weight_sum);
		//							rightvertWtList.clear();
								/*	if(weight_sum==1 && (leftvertWtList.empty()!=true || rightvertWtList.empty()!=true))
									{
										fout<<" e ntv1 "<<ntv1<<" ntv2 "<<ntv2<<endl;	
									}*/
									preSum(leftvertWtList, rightvertWtList, weight_sum, 0, 0, nteWT, postProcout, saveFirst_u1_node, v1_node, ear_premble_vec, 
									mark_eargraph_nodes ,mark_notfree_nodes);
									leftvertWtList.clear();
									rightvertWtList.clear();
									nteWT=0;	

									saveFirst_u1_node=MAX_UINT;	
									weight_sum=0;
								}						

							//printf("b2 %d %d m= %d %d master %d %d\n",store_u1_ear,store_v1_ear,lca[l1].u1,lca[l1].v1,master[j],i);
						}

					}
				}
				if(parent_v1==lca_node)
					break;
				v1_node=parent_v1;
			}
		}
	}		
}
leftvertWtList.clear();
rightvertWtList.clear();
leftvertWtList1.clear();
//#pragma omp end parallel
//==============Final print Master ENDS===========// A001 	A001 

//==============Remove Parallel Edges============
fout.close();
postProcout.close();


//postProcout.open("PostProcessOp.txt",ios::in);
	//std::cout<<" ear_premble_vec size "<<ear_premble_vec.size()<<std::endl;
	//cout<<" Test "<<ear_premble_vec[0]->lnode<<" "<<ear_premble_vec[0]->rnode<<endl;	
	delete [] to_visit;
	delete [] to_visit_threads;
	delete [] visited;
	delete [] thread_backs;
}
