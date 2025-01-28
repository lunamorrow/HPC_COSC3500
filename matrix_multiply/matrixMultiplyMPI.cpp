#include "matrixMultiplyMPI.h"

using namespace std;

void matrixMultiply_MPI(int N, const floatType* A, const floatType* B, floatType* C, int* flags, int flagCount)
{
    MPI_Init(&flagCount, &flags)

    memset(C, 0, sizeof(floatType)*N*N);
    const int block = 256; // best size out of 64, 128, 256 and 512

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int start = my_rank * N/world_size;
    int end = (my_rank+1) * N/world_size;

    __m256d a1, b1, sum1, mul1, a2, b2, sum2, mul2, a3, b3, sum3, mul3, a4, b4, sum4, mul4;
    __m128d b_1, b_2, b_3, b_4;

    int i, k, j, kk, jj, ii, K, J, I1, I2, I3, I4;

    #pragma omp parallel private(a1, b1, sum1, mul1, a2, b2, sum2, mul2, a3, b3, sum3, mul3, a4, b4, sum4, mul4, b_1, b_2, b_3, b_4, i, j, k, jj, kk, ii, J, K, I1, I2, I3, I4) shared(A, B, C)
    {
        #pragma omp for collapse(2)
        for (j = start; j < end; j += block)
        {
            for (k = 0; k < N; k += block)
            {
                for (i = 0; i < N; i += block)
                {
                    for (kk = 0; kk < block; kk++)
                    {
                        K = k + kk;
                        for (jj = 0; jj < block; jj++)
                        {
                            J = j + jj;
                            for (ii = 6; ii < block; ii+=8)
                            {
                                // SECTION 1 //
                                I1 = i + ii - 6;
                                a1 = _mm256_load_pd((double*)&A[(K*N) + I1]);
                                b_1 = _mm_load_pd((double*)&B[(J*N) + K]);
                                b1 = _mm256_broadcast_pd(&b_1);

                                mul1 = _mm256_mul_pd(b1, a1);
                                sum1 = _mm256_load_pd((double*) &C[(J*N) + I1]);
                                sum1 = _mm256_add_pd(sum1, mul1);
                                _mm256_store_pd((double*) &C[(J*N) + I1], sum1);


                                // SECTION 2 //
                                I2 = i + ii - 4;
                                a2 = _mm256_load_pd((double*)&A[(K*N) + I2]);
                                b_2 = _mm_load_pd((double*)&B[(J*N) + K]);
                                b2 = _mm256_broadcast_pd(&b_2);

                                mul2 = _mm256_mul_pd(b1, a2);
                                sum2 = _mm256_load_pd((double*) &C[(J*N) + I2]);
                                sum2 = _mm256_add_pd(sum2, mul2);
                                _mm256_store_pd((double*) &C[(J*N) + I2], sum2);

                                // SECTION 3 //
                                I3 = i + ii - 2;
                                a3 = _mm256_load_pd((double*)&A[(K*N) + I3]);
                                b_3 = _mm_load_pd((double*)&B[(J*N) + K]);
                                b3 = _mm256_broadcast_pd(&b_3);

                                mul3 = _mm256_mul_pd(b1, a3);
                                sum3 = _mm256_load_pd((double*) &C[(J*N) + I3]);
                                sum3 = _mm256_add_pd(sum3, mul3);
                                _mm256_store_pd((double*) &C[(J*N) + I3], sum3);

                                // SECTION 4 //
                                I4 = i + ii;
                                a4 = _mm256_load_pd((double*)&A[(K*N) + I4]);
                                b_4 = _mm_load_pd((double*)&B[(J*N) + K]);
                                b4 = _mm256_broadcast_pd(&b_4);

                                mul4 = _mm256_mul_pd(b1, a4);
                                sum4 = _mm256_load_pd((double*) &C[(J*N) + I4]);
                                sum4 = _mm256_add_pd(sum4, mul4);
                                _mm256_store_pd((double*) &C[(J*N) + I4], sum4);
                            }
                        }
                    }
                }
            }
        }
    }

    MPI_Request request;
    if (my_rank == 0){
        MPI_Irecv(&C[2097152], 4194304, MPI_DOUBLE, 1,
              0, MPI_COMM_WORLD, &request);
	MPI_Send(&C[0], 4194304, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	MPI_Wait(&request, MPI_STATUS_IGNORE);
    } else {
        MPI_Send(&C[2097152], 4194304, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	MPI_Irecv(&C[0], 4194304, MPI_DOUBLE, 0,
              0, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, MPI_STATUS_IGNORE);
    }

    MPI_Finalize()
}
