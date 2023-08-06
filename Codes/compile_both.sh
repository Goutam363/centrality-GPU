nvcc -w -o apsp_both -O3 -arch sm_21 -Xcompiler -fopenmp gpu_apsp.cu graph.cpp apsp_graph.cpp bicc.cpp main.cpp pendant_graph.cpp modified_apsp.cpp hybrid_apsp.cpp
