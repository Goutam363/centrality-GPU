using namespace std;
//#include <fstream>
#include <set>
#include <vector>
//#include <>
//#include <queue>
//#include <map>
//#include "parse.h"
#include "ear_graph.h"


void print(int *srcs,int *dsts,ear_premble**& ear_info,int n,int m){

	ear_premble *obj;
	int k=0,size;
	for(int i=0;i<n-1;i++){

		if(i==n-1)
			size = m;
		else	
			size = srcs[i+1];

		
		for(int j=srcs[i];j<size;j++){
	
			obj = ear_info[j];
			if(obj==NULL)
				continue;
			cout<<i<<" "<<dsts[j]<<endl;
			cout<<obj->lnode <<" "<<obj->rnode<<" "<<" "<<obj->cnt<<" "<<obj->td<<endl;
			k=0;
			while(k<obj->cnt){

				cout<<obj->ld[k]<<" "<<obj->nodes[k]<<" "<<obj->rd[k]<<endl;
				k++;
			}
		}	
	}


}

void show(int *srcs,int* dsts,int n,int m){

	for(int i=0;i<n;i++)
      printf("out_degree_list %d = %d \n",i,srcs[i]);
    for(int i=0;i<m;i++)
      printf("out_array %d = %d \n",i,dsts[i]);
}
void fill(int lnode,int rnode,int *srcs,int *dsts,ear_premble*& obj,ear_premble**&ear_info,int n,int m){

	int size,i;
	int u = lnode;
	int v = rnode;
	for(i=srcs[u];i<srcs[u+1];i++){

		if(dsts[i]==v)
		{
			if(ear_info[i]==NULL)
			{
				ear_info[i] = obj;
				
			}else
			{
				int *nld = new int[ear_info[i]->cnt + obj->cnt];
				int *nrd = new int[ear_info[i]->cnt + obj->cnt];
				int *n_nodes = new int[ear_info[i]->cnt + obj->cnt];

				int k=0;
				for(int j=0;j<ear_info[i]->cnt;j++)
				{
					nld[k] = ear_info[i]->ld[j];
					nrd[k] = ear_info[i]->rd[j];
					n_nodes[k] = ear_info[i]->nodes[j];
					k++;
				}
				for(int j=0;j<obj->cnt;j++)
				{
					nld[k] = obj->ld[j];
					nrd[k] = obj->rd[j];
					n_nodes[k] = obj->nodes[j];
					k++;	
				}
				ear_info[i]->cnt += obj->cnt;
				delete [] ear_info[i]->ld;
				delete [] ear_info[i]->rd;
				delete [] ear_info[i]->nodes;

				ear_info[i]->ld = nld;
				ear_info[i]->rd = nrd;
				ear_info[i]->nodes = n_nodes;
				/*
				ear_info[i]->td is not updated in such cases because td will not be same. This doesnt affect further program as we aren't using td.
				*/
			}
			break;
		}
	}

}

void rearrange(std::vector<int>&req_nodes,std::set<int>&ear_vertices,int n)
{
	std::cout<<"Inside rearrange "<<ear_vertices.size()<<std::endl;
	for(int i=0;i<n;i++)
	{
		//printf("req_nodes = %d\n",i);
		if(ear_vertices.find(i)==ear_vertices.end())
		{
			req_nodes.push_back(i);
			
		}	
	}
	std::cout<<"Required vertices count "<<req_nodes.size()<<std::endl;
}

void read_ear_info(vector<ear_premble *>ear_premble_vec, int* srcs, int* dsts, int n, int m, ear_premble**& ear_info, ear_premble**& self_ear,
std::vector<int>&disc_time,int &active_vert_count)
{

	string line;	
	int flag,lnode,rnode,cnt,td,i=0,ld,rd,nodes,index=0;

	std::set<int> onlyearends;
	//std::set<int> deg2nodes;

	ear_info = new ear_premble*[2*m];
	self_ear = new ear_premble*[n];

	for(i=0;i<(2*m);i++)
	{
		ear_info[i] = NULL;
	}

	for(i=0;i<n;i++)
	{
		self_ear[i] = NULL;
	}
	
	while(index < ear_premble_vec.size())
	{
		
		lnode = ear_premble_vec[index]->lnode;
		rnode = ear_premble_vec[index]->rnode;
		cnt   = ear_premble_vec[index]->cnt;
		td    = ear_premble_vec[index]->td;

		if(cnt>0)
		{	
			onlyearends.insert(lnode);	
			onlyearends.insert(rnode);	
		}	
		i=0;
		if(!(disc_time[lnode]<disc_time[rnode]))
		{
			int tmp = lnode;
			ear_premble_vec[index]->lnode = rnode;
			ear_premble_vec[index]->rnode = tmp;

			while(cnt>0){
		
				tmp = ear_premble_vec[index]->rd[i];
				ear_premble_vec[index]->rd[i] = ear_premble_vec[index]->ld[i];
				ear_premble_vec[index]->ld[i] = tmp;
//				printf("deg2 node = %d ", ear_premble_vec[index]->nodes[i]);
				i++;
				cnt--;
			}
		}
//		printf("\n\n");

		if(lnode == rnode )
		{
			self_ear[lnode] = ear_premble_vec[index];	
		}	
		else
			fill(ear_premble_vec[index]->lnode,ear_premble_vec[index]->rnode,srcs,dsts,ear_premble_vec[index],ear_info,n,(2*m));
		//ear_info[index] = obj;
		
		index++;

	}


	//cout<<" Total vertices in ear graph "<<n<<endl;
	//cout<<" Total vertice degree == 2 in original graph found from ear graph are "<<deg2nodes.size()<<endl;
	//cout<<" Total vertice degree > 2 in ear graphs are "<<onlyearends.size()<<endl;
	active_vert_count += onlyearends.size();
	onlyearends.clear();
//	std::cout<<" end of read_ear_info "<<std::endl;

}
