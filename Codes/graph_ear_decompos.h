#include <cstdio>
#include <vector>
#include <iostream>
#include <cstdlib>

class graph_ear_decompos{
public:
  graph_ear_decompos() 
  {
    num_verts = 0;
    num_edgs = 0;
  }

  ~graph_ear_decompos() {};

void check_connect(int n, int m, unsigned int *&R, int *&C, bool *&visited,int s)
{
    int end;
    visited[s] = true;
    if(s==n-1)
      end = m;
    else
      end = R[s+1];
    for(int f = R[s]; f < end; f++)
    {
        if(visited[C[f]]==false)
         {
          check_connect(n,m,R,C,visited,C[f]);
         } 
    } 


}

  // create csr from edgelist 
void init(int n, unsigned int m, int* srcs, int* dsts,int *wgt)
{
    num_verts = n;
    num_edgs = m;   
    avg_out_degree = 0;
    
    out_array = new int[num_edgs];
    weights = new int[num_edgs];        //A001
    out_degree_list = new unsigned int[num_verts+1];
    out_degree_list[0] = 0;
    non_tree_edge=new bool[num_edgs]; //A001
    serial=new unsigned int[num_edgs];         //A001
    master=new unsigned int[num_edgs];         //A001
    lca_array=new unsigned int[num_edgs];
    bridges=new unsigned int[num_edgs];         //A001
    bool *visited = new bool[num_verts];

    // get out verts
    unsigned int* temp_counts = new unsigned int[num_verts];
    unsigned int* degree_array= new unsigned int[num_verts]; 	//A001	Stores degree of every vertex
    unsigned int var1;
    unsigned int degree2Vertices=0,degree3Vertices=0,degree1Vertices=0,degree0Vertices=0;
    for (int i = 0; i < num_verts; ++i)
    {  
      temp_counts[i] = 0;
      visited[i] = false;
    } 
    for (unsigned int i = 0; i < num_edgs; ++i)
    {
      ++temp_counts[srcs[i]];
      non_tree_edge[i]=0;         //initialized to 0 A001
      lca_array[i]=999999999;              //initialized to 0 A001
      serial[i]=i;                 //initialized to 0 A001
      master[i]=999999999;                 //initialized to 0 A001
    }
  //  printf("hi2\n");   //A001
    for (int i = 0; i < num_verts; ++i)
      if(temp_counts[i]!=0)
        out_degree_list[i+1] = out_degree_list[i] + temp_counts[i];

    for(int i=0;i<num_verts;i++)	//A001		Copy degrees to degree array
    	degree_array[i]=temp_counts[i];       //A001 
    
	  std::copy(out_degree_list, out_degree_list + num_verts, temp_counts);
    for (unsigned int i = 0; i < num_edgs; ++i)
    {
      var1=temp_counts[srcs[i]]++;
      out_array[var1] = dsts[i];
      weights[var1] = wgt[i];     //A001
    }
    delete [] temp_counts;

    long unsigned degree_sum = 0;
    max_out_degree = 0;
    max_degree_vert = 0;
    for (int i = 0; i < num_verts; ++i)
    {
      int degree = out_degree_list[i+1] - out_degree_list[i];
	    degree_sum += (long unsigned)degree;
      if (degree > max_out_degree)
      {
        max_out_degree = degree;
        max_degree_vert = i;
      }
    }
    avg_out_degree = (double)degree_sum / (double)num_verts;

#if GET_GRAPH_PARAMS
    printf("n: %d\n", n);
    printf("m: %d\n", m / 2);
    printf("avg: %9.6f\n", avg_out_degree);
    printf("max: %d\n", max_out_degree);
    //printf("dia: %d\n", calc_diameter());
#endif
   // printf("max_degree_vert=%d\n",max_degree_vert); // 30-9   A001 
	//To print degrees of vertices
	//To count 2 degree and degree > 3 vertices
	for(int j=0;j<num_verts;j++)
	{
		if(degree_array[j]==0)
			degree0Vertices++;		
		if(degree_array[j]==1)
			degree1Vertices++;		
		if(degree_array[j]==2)
			degree2Vertices++;
		if(degree_array[j]>2)
			degree3Vertices++;		
		//	cout<<" "<<degree_array[j]<<" ";			
	}
	//printf(" Degree 0 count is %u\n",degree0Vertices);
 	//printf(" Degree 1 count is %u\n",degree1Vertices);
	//printf(" Degree 2 count is %u\n",degree2Vertices);
	//printf(" Degree greater than 2 count is %u\n",degree3Vertices);//" "<<max_count<<" "<<pendantcnt<<endl;
  //   printf("Degree 2 n >2 Vertices d(1) %d %d %d %d\n",degree2Vertices,degree3Vertices,degree1Vertices,degree0Vertices); 
   //exit(0);  
#if CDEBUG
    //===Print CSR format out_array indicates destination======
    for (unsigned int i = 0; i <num_edgs; ++i)
        printf(" out_array=%d %d %d\n",out_array[i],wgt[i],serial[i]);
    //===========out_degree_list indiactes source or starting position=====
    for (int i = 0; i < num_verts; ++i)
         printf(" o_l=%d \n",out_degree_list[i]);
#endif


      check_connect(n,m,out_degree_list,out_array, visited, 0);
      for(int f=0; f < n; f++)
       {
          if(visited[f]==false)
          {
            printf("Graph is disconnected \n");
            exit(0);
          }
       } 
}

  /*int calc_diameter()
  {
    int* queue = new int[num_verts];
    int front = 0;
    int back = 0;
    int max_dist = 0;
    int abs_max_dist = 0;

    bool* visited = new bool[num_verts];
    for (int i = 0; i < num_verts; ++i)
      visited[i] = false;

    srand(time(0));

    int cur_dist;
    for (int q = 0; q < 5; ++q)
    {
      cur_dist = 0;
      max_dist = 0;

      queue[0] = (int)(((double) num_verts) * 
              ((double) rand()/(((double) RAND_MAX) + 1.0)));
      visited[queue[0]] = true;
      front = 0;
      back = 1;
      bool go = true;

      while (go)
      {
        while (back - front)
        {
          int size = back - front;
          for (int i = 0; i < size; ++i)
          {
            int vert = queue[front + i];
            int out_deg = out_degree(vert);
            int* outs = out_vertices(vert);
            int out;
            for (int i = 0; i < out_deg; ++i)
            {
              out = outs[i];
              if (!visited[out])
              {
                visited[out] = true;
                queue[back++] = out;
              }
            }
          }
          ++cur_dist;
          front += size;
        }
        --front;
        for (int i = 0; i < num_verts; ++i)
          visited[i] = false;

        if (cur_dist > max_dist)
        {
          max_dist = cur_dist;
          cur_dist = 0;
          go = true;
          queue[0] = queue[front];
          front = 0;
          back = 1;
          visited[queue[0]] = true;
        }
        else
          go = false;
      }

      if (max_dist > abs_max_dist)
        abs_max_dist = max_dist;
    }

    delete [] queue;
    delete [] visited;

    return abs_max_dist;
  }*/
  
  double get_avg_out_degree()
  {
    return avg_out_degree;
  }

  int get_max_out_degree()
  {
    return max_out_degree;
  }

  int get_max_degree_vert()
  {
    return max_degree_vert;
  }

  int* get_out_vertices()
  {
    return out_array;
  }
  int* get_out_weights()   //A001
  { 
    return weights;       //A001
  }

  bool* get_out_non_tree_edge() //A001
  {
    return non_tree_edge;       //A001
  }

  unsigned int* get_master() //A001
  {
    return master;       //A001
  }
  unsigned int* get_serial() //A001
  {
    return serial;       //A001
  }
  unsigned int* get_lca_array() //A001
  {
    return lca_array;       //A001
  }

  int* out_vertices(int v)
  {
    return (&out_array[out_degree_list[v]]);
  }

  int out_degree(int v)
  {
    return out_degree_list[v+1] - out_degree_list[v];
  }

  int* outs() const
  {
    return out_array;
  }



  unsigned int* get_out_degree_list() const
  {
    return out_degree_list;
  }

  int num_vertices() const
  {
    return num_verts;
  }
  
  unsigned int num_edges() const
  {
    return num_edgs;
  }

  graph_ear_decompos& operator= (const graph_ear_decompos& param)
  {
    if (num_edgs)
      delete [] out_array;
    if (num_verts)
      delete [] out_degree_list;

    num_verts = param.num_vertices();
    num_edgs = param.num_edges();    
    out_array = new int[num_edgs];
    out_degree_list = new unsigned int[num_verts+1]; 
    avg_out_degree = param.avg_out_degree;   
    std::copy(param.outs(), param.outs() + num_edgs, out_array);
    std::copy(param.get_out_degree_list(), param.get_out_degree_list() + (num_verts+1), out_degree_list);
      
    return *this;
  }
  
  void clear()
  {
    if (num_edgs)
      delete [] out_array;
    if (num_verts)
      delete [] out_degree_list;

    num_verts = 0;
    num_edgs = 0;
  }
  
private:
  int num_verts;
  unsigned int num_edgs;
  double avg_out_degree;
  int max_out_degree;
  int max_degree_vert;
  int* out_array;
  unsigned int* out_degree_list;
  bool *non_tree_edge; //A001
  unsigned int *serial; //A001
  unsigned int *master; //A001
  unsigned int *lca_array;
  unsigned int *bridges; //A001
  int* weights;     //A001
};
