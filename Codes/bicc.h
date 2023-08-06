/** 
 *
 * Class for finding the biconnected components of the given input graph, as edge lists;
 *	
 *	Init 	: 	BiCC bicc = BiCC(filename/Graph);
 *	Usage 	: 	vector<edgelist> edgelists = bicc.run_bicc();
 *			Returns each BiCC as an edgelist;
 **/

#ifndef BICC_H_
#define BICC_H_

#include "graph.h"
#include <stack>
#include <map>
#include <vector>
#include <utility>
//map<int,map<int,int> > bccvertno;					//Block new vertex old vertex mapping
//map<int,vii> bccvertno;
//vector<int> blocklevel;						//To store BCC levels
class BiCC {
	private:
		//int prv_u;	
		//int bcc_cnt;
	public :
		Graph graph;
		vi dfs_number;
		vi high_water;
		stack<iii> edge_stack;
		//map<int,vviii> bccmap;				//CH01
		map<int,vi> bccmap;					// used to generate block tree
		//vector<bool> art
		//vvi bccvector;
		vviii bccmatrix;					//Matrix to store APSP bet all BCC's
		vi artdegree;						//Stores degree of vertices on block tree	
		map<int,int> artmap;					//Holds block tree art point mapping to its new value
		map<int,int> inv_art;
		map<int,viii> bccrelabel;				//Not used
		map<int,vviii>::iterator bccitr; 			//CH01
		static map<int,map<int,pair<int,int> > > bccvertno;	//Block new vertex old vertex mapping <bccnumber,<relabeld val,pair<oldvertex,reachval>>>
		static map<int,map<int,int> > rev_bccvertno;		//Block oldartpt to new_art mapping <bccnumber,<relabeldval,oldvertex>>
		map <int,int> bccvrtcnt;		//Mapping of BCC number to number of vertices in that BCC
		map<int,vi> bcclevelmap;			//level to BCC's on that level mapping
		static vector<int> blocklevel;				//To store BCC levels
		
		vector<bool> leafnodes;					//To store leaf nodes of block tree		
		vii block_index;					//Stores block id and new variable 
		vii bccnontree;	 					//Store non tree edges in bcc tree
		vi bccparent;						//Holds block tree structure
		vi bccmapping;						//Holds mapping of BCC's to integer value
		vviii subgraphs_edgelist;
		//static vi reach;					//To hold number of nodes reachable thr i'th
		int n_vertices;
		int vertices_count;					//Total number of vertices in block tree
		int last_dfs_number;
		//int n_vertices;
		int prv_u;	
		int bcc_cnt;	
		BiCC();
		BiCC(char *file);
		BiCC(Graph &graph);
		vviii run_bicc();
		void print_bicc_hist_to_file(char *file);
		void bicc_tarjan(int u, int parent );
		void print_edgelist();
		void fillbccparent(int totalvertex);
		void update_op(int art,int n1,int n2);
		void pc1(vii bcc_sizes);
		void printoutput(float *&bc_result);
		void postproc1();
		void process_ij(vi &i_lca,int i,int j,int mbi,int mbj);
		vi reverse_vect(vi i_lca,vi j_lca,int ilcasize,int jlcasize);
		void reachfunction(int totalvertex, float *&bc_result);
		void to_eargraph(int g_nodes,int g_edges, vector<vector<pair<pair<int,int>,int > > >&edgelists, map<int,map<int,pair<int,int> > >&bccvertno 
		,map<int,int>& artmap, int num_bcc, map<int,int>&bccvrtcnt, map<int,int>& inv_artmap,  float *&bc_result, double &g_free_time,
		int &active_vert_count, int &free_vert_count);
};	
#endif
