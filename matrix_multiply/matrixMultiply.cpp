#include "matrixMultiply.h"

void matrixMultiply(int N, const floatType* A, const floatType* B, floatType* C, int* args, int argCount)
{
    memset(C, 0, sizeof(floatType)*N*N);
    const int block = 256; // best out of 64, 128, 256 and 512

    __m256d a1, b1, sum1, mul1, a2, b2, sum2, mul2, a3, b3, sum3, mul3, a4, b4, sum4, mul4; 
    __m128d b_1, b_2, b_3, b_4;

    int i, k, j, kk, jj, ii, K, J, I1, I2, I3, I4;

    #pragma omp parallel private(a1, b1, sum1, mul1, a2, b2, sum2, mul2, a3, b3, sum3, mul3, a4, b4, sum4, mul4, b_1, b_2, b_3, b_4, i, j, k, jj, kk, ii, J, K, I1, I2, I3, I4) shared(A, B, C)
    {
        #pragma omp for collapse(2)
        // 2 gives grade of 5.46-6.07
        // 3 gives grade of 5.50-5.71
        // 4 gives grade of 5.45-5.66
        for (i = 0; i < N; i += block)
        {
            for (k = 0; k < N; k += block)
            {
                for (j = 0; j < N; j += block)
                {
                    for (kk = 0; kk < block; kk++) //try ikj, jki - might be faster (or ikj, kji for c accumulator)
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
                                // printf("a1 indexes = %d and %d\n", (K*N) + I1, (K*N) + I1 + 1);
                                // printf("b1 index = %d\n", (J*N) + K);
                                b_1 = _mm_load_pd((double*)&B[(J*N) + K]);
                                b1 = _mm256_broadcast_pd(&b_1);

                                mul1 = _mm256_mul_pd(b1, a1);
                                // printf("c1 indexes = %d\n", (J*N) + I1, (J*N) + I1 + 1);
                                sum1 = _mm256_load_pd((double*) &C[(J*N) + I1]);
                                sum1 = _mm256_add_pd(sum1, mul1);
                                _mm256_store_pd((double*) &C[(J*N) + I1], sum1);

                                
                                // SECTION 2 //
                                I2 = i + ii - 4;
                                a2 = _mm256_load_pd((double*)&A[(K*N) + I2]);
                                // printf("a2 indexes = %d and %d\n", (K*N) + I2, (K*N) + I2 + 1);
                                // printf("b2 index = %d\n", (J*N) + K);
                                b_2 = _mm_load_pd((double*)&B[(J*N) + K]);
                                b2 = _mm256_broadcast_pd(&b_2);

                                mul2 = _mm256_mul_pd(b1, a2);
                                //  printf("c2 indexes = %d\n", (J*N) + I2, (J*N) + I2 + 1);
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
}
