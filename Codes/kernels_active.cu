#include <cuda.h>
#include <cuda_runtime.h>

// Assuming ear_premble is defined
struct ear_premble {
    int start_node;
    int end_node;
    // Note: std::vector<int> nodes can't be used directly in CUDA, so assume nodes are handled differently
    int* nodes; // Device pointer to node list
    int node_count;
};

// Betweenness centrality kernel for active nodes
__global__ void bc_kernel_active(int* R, int* C, int* F, int g_n, int g_m, float* bc, int* sampled_nodes, int sample_size, ear_premble* ear_data) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= sample_size) return;

    int s = sampled_nodes[idx]; // Source node

    // Shared memory for BFS
    extern __shared__ char sh_mem[];
    int* Q = (int*)sh_mem; // Queue
    int* S = &Q[g_n]; // Stack
    float* delta = (float*)&S[g_n]; // Dependency
    int* sigma = (int*)&delta[g_n]; // Shortest path counts
    int* d = &sigma[g_n]; // Distances

    // Initialize arrays
    for (int i = threadIdx.x; i < g_n; i += blockDim.x) {
        delta[i] = 0.0f;
        sigma[i] = 0;
        d[i] = -1;
    }
    __syncthreads();

    // BFS initialization
    if (threadIdx.x == 0) {
        Q[0] = s;
        d[s] = 0;
        sigma[s] = 1;
    }

    int q_front = 0, q_rear = 1, s_top = 0;
    __syncthreads();

    // BFS
    while (q_front < q_rear) {
        int v = Q[q_front++];
        int start = F[v];
        int end = F[v + 1];
        for (int i = start + threadIdx.x; i < end; i += blockDim.x) {
            int w = R[i];
            if (atomicCAS(&d[w], -1, d[v] + 1) == -1) {
                Q[q_rear++] = w;
                sigma[w] = 0;
            }
            if (d[w] == d[v] + 1) {
                atomicAdd(&sigma[w], sigma[v]);
            }
        }
        __syncthreads();
    }

    // Dependency accumulation (simplified, adjust for ear_data)
    for (int i = threadIdx.x; i < g_n; i += blockDim.x) {
        if (d[i] >= 0 && i != s) S[s_top++] = i;
    }
    __syncthreads();

    while (s_top > 0) {
        int w = S[--s_top];
        for (int i = F[w] + threadIdx.x; i < F[w + 1]; i += blockDim.x) {
            int v = R[i];
            if (d[v] == d[w] - 1) {
                float val = (sigma[v] / (float)sigma[w]) * (1.0f + delta[w]);
                atomicAdd(&delta[v], val);
            }
        }
        if (threadIdx.x == 0 && w != s) {
            atomicAdd(&bc[w], delta[w]);
        }
        __syncthreads();
    }

    // Use ear_data (example: adjust centrality for ear nodes)
    // Note: This is a placeholder; adjust based on actual ear_premble usage
    if (threadIdx.x == 0 && ear_data != nullptr) {
        for (int i = 0; i < sample_size; ++i) {
            if (ear_data[i].start_node == s || ear_data[i].end_node == s) {
                bc[s] += delta[s]; // Example adjustment
            }
        }
    }
}