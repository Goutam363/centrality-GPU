//#include "ear.h"
#include <cstdio>
#include <vector>
#include <utility>
#include <iostream>
#include <queue>
#include <algorithm>
#define THRESHOLD 2000
using namespace std;

void dobfs(int*R,int* F,int* C,int n,int m,vector<int>&parents,vector<int>&levels,vector<int>&disc_time, int min_degree_node)
{
	int s = min_degree_node;
	int start,end;
	vector<bool> mark(n,0);
	queue<int> q;
	q.push(s);
	mark[s] = 1;

	levels[s] = 0;
	//printf("In BFS n = %d\n",n);
	int t=1;
	disc_time[s] = t++;
	while(!q.empty())
	{
		s = q.front();
		q.pop();
		start = R[s]; end = R[s+1];
		//cout<<"start "<<start<<" end "<<end<<endl;	
		for(int i=start;i<end;i++)
		{
			//cout<<C[i]<<" "<<mark[C[i]]<<endl;	
			if(mark[C[i]]==0)
			{
				mark[C[i]] = 1;
				disc_time[C[i]] = t++;
				parents[C[i]] = s;
				levels[C[i]] = (levels[s] + 1);
				//cout<<" lev of "<<s<<" "<<levels[s]<<" lev of "<<C[i]<<" "<<levels[C[i]]<<endl;
				q.push(C[i]);
			}
		}
	}

	mark.clear();
	//q.clear();

	/*for(int i=0;i<levels.size();i++)
	{
		cout<<" lev "<<i<<" = "<<levels[i]<<std::endl;
	}	
	printf("BFS ends\n");*/
}


void arrange_level_order(vector<int> &level_offset,vector<int>&level_order,vector<int>&levels,int n,int max_earlevel,vector<int>&disc_time)
{

	vector<pair<int,int> > vii;
	vector<int>n_level_offset;
	vector<int>n_level_order;
	int size=0;

	//cout<<" arrange_level_order begins "<<endl;
	for(int i=0;i<n;i++)
	{
		if(levels[i]!=-1)
		{
			vii.push_back(make_pair(levels[i],i));		
			//cout<<" i = "<<i<<" "<<levels[i]<<" disc_time "<<disc_time[i]<<endl;
		}	
	}

	sort(vii.begin(),vii.end());

	level_offset[0] = 0;
	for(int i=0;i<vii.size();i++)
	{
		level_offset[vii[i].first+1]++;
		level_order[i] = vii[i].second;

	}	

	size = level_offset.size();
	int tmp=0;
	n_level_offset.push_back(0);
	for(int i=1;i<size;i++)
	{

		if((tmp + level_offset[i])<THRESHOLD)
			tmp += level_offset[i];
		else
		{
			n_level_offset.push_back(tmp);
			tmp=level_offset[i];
		}
	}

	n_level_offset.push_back(tmp);

	size = n_level_offset.size();
	for(int i=1;i<size;i++)
	{
		n_level_offset[i] += n_level_offset[i-1]; 

	}

	level_offset.clear();
	level_offset = n_level_offset;
	
	size = level_offset.size();
	vii.clear();
	/*cout<<" level_offset "<<endl;
	for(int i=0;i<size;i++)
	{
		cout<<" i = "<<level_offset[i]<<endl;
	}

	cout<<endl;
	
	int k;
	for(int i=0;i<size;i++)
	{
		cout<<" level is "<<i<<endl;
		if(i==size-1)
			k= size;
		else
			k = level_offset[i+1];
		for(int j=level_offset[i];j<k;j++)
		{
			cout<<level_order[j]<<" ";
		}
		cout<<endl;
	}


	cout<<"arrange_level_order ends"<<endl;	
	
	cout<<" n_level_offset is "<<endl;
	for(int i=0;i<n_level_offset.size();i++)
	{
		cout<<"level is "<<i<<" nodes are "<<n_level_offset[i]<<endl;
	}*/
	return;

}
