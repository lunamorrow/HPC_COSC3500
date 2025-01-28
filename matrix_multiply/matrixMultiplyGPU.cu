#include "matrixMultiplyGPU.cuh"

__host__ void matrixMultiply_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C, int* flags, int flagCount)
{

    // this is the CPU, it hosts the GPU device -> only calls CPU functions
    // call the Kermal function!
    // declare the number of blocks per grid and the number of threads per block
    // int num_b = N/32; //((N+3)/4 + 32 - 1);
    dim3 num_blocks = dim3(64, 64, 1);
    dim3 num_threads = dim3(32, 32, 1);
    
    const size_t sz = N * N * sizeof(floatTypeCUDA);
    cudaMemset(C, 0, sz);
    
    matrixMultiplyKernel_GPU<<<num_blocks, num_threads, sizeof(floatTypeCUDA)*1024>>>(N, A, B, C, 0, 0, 0);
}


__global__ void matrixMultiplyKernel_GPU(int N, const floatTypeCUDA* A, const floatTypeCUDA* B, floatTypeCUDA* C, int flag0, int flag1, int flag2)
{
    int tid = threadIdx.x * blockDim.y + threadIdx.y;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    extern __shared__ double smem[];
    
    smem[tid] = __uint2double_rn(0);
        for (int k = 0; k < N; k+=16){
	    smem[tid] += ((B[col * N + k].x * A[k * N + row].x) +
	            (B[col * N + k + 1].x * A[(k+1) * N + row].x) +
	    (B[col * N + k+2].x * A[(k+2) * N + row].x) +
	    (B[col * N + k+3].x * A[(k+3) * N + row].x) +
	    (B[col * N + k+4].x * A[(k+4) * N + row].x) +
	    (B[col * N + k+5].x * A[(k+5) * N + row].x) +
	    (B[col * N + k+6].x * A[(k+6) * N + row].x) +
	    (B[col * N + k+7].x * A[(k+7) * N + row].x) +
	    (B[col * N + k+8].x * A[(k+8) * N + row].x) +
	    (B[col * N + k+9].x * A[(k+9) * N + row].x) +
	    (B[col * N + k+10].x * A[(k+10) * N + row].x) +
	    (B[col * N + k+11].x * A[(k+11) * N + row].x) +
	    (B[col * N + k+12].x * A[(k+12) * N + row].x) +
	    (B[col * N + k+13].x * A[(k+13) * N + row].x) +
	    (B[col * N + k+14].x * A[(k+14) * N + row].x) +
	    (B[col * N + k+15].x * A[(k+15) * N + row].x));


        }
    
    C[col * N + row].x = smem[tid];

}
