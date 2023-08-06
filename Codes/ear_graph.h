#include "graph_ear_decompos.h"
#include <assert.h>
class ear_premble
{
public:
	
	int lnode;
	int rnode;
	int cnt;
	int td;
	int *ld,*rd,*nodes;

	ear_premble(){}
	ear_premble(int lnode,int rnode,int cnt,int td){

		this->lnode = lnode;
		this->rnode = rnode;
		this->cnt = cnt;
		this->td = td;
		ld = new int[cnt];
		rd = new int[cnt];
		nodes = new int[cnt];


	}
	~ear_premble()
	{
		delete [] ld;
		delete [] rd;
		delete [] nodes;
	}
};

void bcc_bfs(graph_ear_decompos& g, std::vector<ear_premble*>&ear_premble_vec, std::vector<bool>& mark_eargraph_nodes, 
std::vector<bool>& mark_notfree_nodes);

void parse_active_nodes(std::vector<ear_premble*>&ear_premble_vec, int n, int *&R, int *&C, int*&F, int &ear_m, int &min_degree_node);

void arrange_level_order(std::vector<int> &level_offset, std::vector<int>&level_order, std::vector<int>&levels, int n, int max_earlevel, 
std::vector<int>&disc_time);

void read_ear_info(std::vector<ear_premble*>ear_premble_vec, int* srcs, int* dsts, int n, int m, ear_premble**& ear_info, ear_premble**& self_ear,
std::vector<int>&disc_time, int &active_vert_count);

void begin_gpu(int *&R, int *&C, int *&F, int n, int m, ear_premble **&ear_info, ear_premble**&self_ear, std::vector<int>&req_nodes, 
std::vector<int>&artnodes, float *&bc_gpu_free, double &g_free_time);

void begin_gpu(int *&R, int *&C, int *&F, int n, int m, ear_premble **&ear_info, std::vector<int>&parents, std::vector<int>&levels, std::vector<int>&level_order, 
std::vector<int>&level_offset, int max_earlevel, std::vector<int>&artnodes, float *&bc_gpu, double &g_free_time, int *&ear_R, int *&ear_C);

void dobfs(int*R,int* F,int* C,int n,int m,std::vector<int>&parents,std::vector<int>&levels,std::vector<int>&disc_time, int min_degree_node);