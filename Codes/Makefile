#CC=icpc
NCC=nvcc

all:
       
        #g++ bcc.cpp -I ../include/ -fopenmp -O3 -o bcc 
		$(NCC) -w -g -lcudadevrt -lineinfo -o betc  -arch=sm_35 -Xcompiler -fopenmp begin_gpupgm.cu kernels_free.cu kernels_active.cu postprocesses.cu graph.cpp bicc.cpp main.cpp pendant_graph.cpp modified_apsp.cpp arrangelevels.cpp bcc_bfs.cpp readeardata.cpp to_eargraph_construct.cpp make_ear_graph.cpp
        
clean:
	rm apsp
