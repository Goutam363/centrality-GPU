#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>

// Assuming ear_premble is defined elsewhere (e.g., in a header)
struct ear_premble {
    int start_node;
    int end_node;
    std::vector<int> nodes;
    // Add other fields as needed
};

// Kernel declarations
__global__ void bc_kernel_active(int* R, int* C, int* F, int g_n, int g_m, float* bc, int* sampled_nodes, int sample_size, ear_premble* ear_data);
__global__ void bc_kernel_free(int* R, int* C, int* F, int g_n, int g_m, float* bc, int* sampled_nodes, int sample_size, ear_premble* ear_data);

// First begin_gpu overload (called at line 165)
void begin_gpu(int*& R, int*& C, int*& F, int g_n, int g_m, ear_premble**& ear_active, ear_premble**& ear_free,
               std::vector<int>& active_nodes, std::vector<int>& free_nodes, float*& h_bc, double& gpu_time) {
    // Start timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    // Determine sample size
    int sample_size = static_cast<int>(std::sqrt(static_cast<float>(g_n)));
    if (sample_size < 1) sample_size = 1;

    // Sample nodes (prefer active nodes, fallback to all nodes)
    std::vector<int> nodes = active_nodes.empty() ? std::vector<int>(g_n) : active_nodes;
    if (nodes.size() == g_n) {
        for (int i = 0; i < g_n; ++i) nodes[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(nodes.begin(), nodes.end(), g);
    if (sample_size > nodes.size()) sample_size = nodes.size();
    std::vector<int> sampled_nodes(nodes.begin(), nodes.begin() + sample_size);

    // Device memory
    int *d_R, *d_C, *d_F, *d_sampled_nodes;
    float* d_bc;
    ear_premble *d_ear_active, *d_ear_free;

    // Allocate device memory
    cudaMalloc(&d_R, g_m * sizeof(int));
    cudaMalloc(&d_C, g_m * sizeof(int));
    cudaMalloc(&d_F, (g_n + 1) * sizeof(int));
    cudaMalloc(&d_bc, g_n * sizeof(float));
    cudaMalloc(&d_sampled_nodes, sample_size * sizeof(int));
    cudaMalloc(&d_ear_active, sizeof(ear_premble) * active_nodes.size());
    cudaMalloc(&d_ear_free, sizeof(ear_premble) * free_nodes.size());

    // Copy data to device
    cudaMemcpy(d_R, R, g_m * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, g_m * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F, F, (g_n + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemset(d_bc, 0, g_n * sizeof(float));
    cudaMemcpy(d_sampled_nodes, sampled_nodes.data(), sample_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ear_active, *ear_active, sizeof(ear_premble) * active_nodes.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ear_free, *ear_free, sizeof(ear_premble) * free_nodes.size(), cudaMemcpyHostToDevice);

    // Kernel launch
    int threadsPerBlock = 256;
    int blocks = (sample_size + threadsPerBlock - 1) / threadsPerBlock;
    bc_kernel_active<<<blocks, threadsPerBlock>>>(d_R, d_C, d_F, g_n, g_m, d_bc, d_sampled_nodes, sample_size, d_ear_active);
    cudaDeviceSynchronize();
    bc_kernel_free<<<blocks, threadsPerBlock>>>(d_R, d_C, d_F, g_n, g_m, d_bc, d_sampled_nodes, sample_size, d_ear_free);
    cudaDeviceSynchronize();

    // Copy results back
    cudaMemcpy(h_bc, d_bc, g_n * sizeof(float), cudaMemcpyDeviceToHost);

    // Scale results
    float scaling_factor = static_cast<float>(g_n) / sample_size;
    for (int i = 0; i < g_n; ++i) {
        h_bc[i] *= scaling_factor;
    }

    // Record time
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    gpu_time = milliseconds / 1000.0;

    // Free memory
    cudaFree(d_R);
    cudaFree(d_C);
    cudaFree(d_F);
    cudaFree(d_bc);
    cudaFree(d_sampled_nodes);
    cudaFree(d_ear_active);
    cudaFree(d_ear_free);
}

// Second begin_gpu overload (called at line 174)
void begin_gpu(int*& R, int*& C, int*& F, int g_n, int g_m, ear_premble**& ear_active, std::vector<int>& active_nodes,
               std::vector<int>& free_nodes, std::vector<int>& ear_active_nodes, std::vector<int>& ear_free_nodes,
               int ear_count, std::vector<int>& ear_nodes, float*& h_bc, double& gpu_time, int*& ear_start, int*& ear_end) {
    // Start timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    // Determine sample size
    int sample_size = static_cast<int>(std::sqrt(static_cast<float>(g_n)));
    if (sample_size < 1) sample_size = 1;

    // Sample nodes
    std::vector<int> nodes = active_nodes.empty() ? std::vector<int>(g_n) : active_nodes;
    if (nodes.size() == g_n) {
        for (int i = 0; i < g_n; ++i) nodes[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(nodes.begin(), nodes.end(), g);
    if (sample_size > nodes.size()) sample_size = nodes.size();
    std::vector<int> sampled_nodes(nodes.begin(), nodes.begin() + sample_size);

    // Device memory
    int *d_R, *d_C, *d_F, *d_sampled_nodes, *d_ear_start, *d_ear_end;
    float* d_bc;
    ear_premble* d_ear_active;

    // Allocate device memory
    cudaMalloc(&d_R, g_m * sizeof(int));
    cudaMalloc(&d_C, g_m * sizeof(int));
    cudaMalloc(&d_F, (g_n + 1) * sizeof(int));
    cudaMalloc(&d_bc, g_n * sizeof(float));
    cudaMalloc(&d_sampled_nodes, sample_size * sizeof(int));
    cudaMalloc(&d_ear_active, sizeof(ear_premble) * ear_count);
    cudaMalloc(&d_ear_start, ear_count * sizeof(int));
    cudaMalloc(&d_ear_end, ear_count * sizeof(int));

    // Copy data to device
    cudaMemcpy(d_R, R, g_m * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C, C, g_m * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_F, F, (g_n + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemset(d_bc, 0, g_n * sizeof(float));
    cudaMemcpy(d_sampled_nodes, sampled_nodes.data(), sample_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ear_active, *ear_active, sizeof(ear_premble) * ear_count, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ear_start, ear_start, ear_count * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ear_end, ear_end, ear_count * sizeof(int), cudaMemcpyHostToDevice);

    // Kernel launch
    int threadsPerBlock = 256;
    int blocks = (sample_size + threadsPerBlock - 1) / threadsPerBlock;
    bc_kernel_active<<<blocks, threadsPerBlock>>>(d_R, d_C, d_F, g_n, g_m, d_bc, d_sampled_nodes, sample_size, d_ear_active);
    cudaDeviceSynchronize();
    // Note: Free nodes may be handled differently based on ear_nodes
    bc_kernel_free<<<blocks, threadsPerBlock>>>(d_R, d_C, d_F, g_n, g_m, d_bc, d_sampled_nodes, sample_size, d_ear_active);
    cudaDeviceSynchronize();

    // Copy results back
    cudaMemcpy(h_bc, d_bc, g_n * sizeof(float), cudaMemcpyDeviceToHost);

    // Scale results
    float scaling_factor = static_cast<float>(g_n) / sample_size;
    for (int i = 0; i < g_n; ++i) {
        h_bc[i] *= scaling_factor;
    }

    // Record time
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    gpu_time = milliseconds / 1000.0;

    // Free memory
    cudaFree(d_R);
    cudaFree(d_C);
    cudaFree(d_F);
    cudaFree(d_bc);
    cudaFree(d_sampled_nodes);
    cudaFree(d_ear_active);
    cudaFree(d_ear_start);
    cudaFree(d_ear_end);
}

// Error checking utility
void checkCudaError(cudaError_t err, const char* msg) {
    if (err != cudaSuccess) {
        std::cerr << msg << ": " << cudaGetErrorString(err) << std::endl;
        exit(1);
    }
}