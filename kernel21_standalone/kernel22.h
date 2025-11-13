#ifndef KERNEL22_H
#define KERNEL22_H

#include <stdlib.h>

// Access macros for matrices in column-major format
#define A(i,j) A[(i)+(j)*LDA]
#define B(i,j) B[(i)+(j)*LDB]
#define C(i,j) C[(i)+(j)*LDC]

#define min(a, b) (((a) < (b)) ? (a) : (b))

// Cache blocking parameters - tuned for AMD EPYC with:
// L1D: 48KB per core
// L2: 1MB per core
// L3: 32MB per CCX
#define M_BLOCK 96   // ((48KB/3)/8bytes) ≈ 96 rows for L1D
#define N_BLOCK 256  // Utilize more of L2 for B blocks
#define K_BLOCK 384  // Larger K allows more reuse in L2

// Register blocking sizes (4x4 remains optimal for scalar code)
#define MR 4
#define NR 4

// Scale matrix C by beta
void scale_c_k22(double *C, int M, int N, int LDC, double beta) {
    if (beta == 0.0) {
        for (int j = 0; j < N; j++) { // not vectorized use -vec-threshold0 to override. unrolled
            for (int i = 0; i < M; i++) {
                C(i,j) = 0.0;
            }
        }
    } else if (beta != 1.0) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) { // was vectorized, length=8
                C(i,j) *= beta;
            }
        }
    }
}

// Optimized micro-kernel for 4x4 blocks
static inline void micro_kernel_k22(int k_size,
                                  double alpha,
                                  const double *A,
                                  const double *B,
                                  double *C,
                                  int LDC) {
    // Load and initialize accumulator registers
    double c00 = C(0,0), c01 = C(0,1), c02 = C(0,2), c03 = C(0,3);
    double c10 = C(1,0), c11 = C(1,1), c12 = C(1,2), c13 = C(1,3);
    double c20 = C(2,0), c21 = C(2,1), c22 = C(2,2), c23 = C(2,3);
    double c30 = C(3,0), c31 = C(3,1), c32 = C(3,2), c33 = C(3,3);

    // Main computation loop
    for (int k = 0; k < k_size; k++) {  // was vectorized length=8
        double a0 = alpha * A[0 + k*MR];
        double a1 = alpha * A[1 + k*MR];
        double a2 = alpha * A[2 + k*MR];
        double a3 = alpha * A[3 + k*MR];

        double b0 = B[k + 0*k_size];
        double b1 = B[k + 1*k_size];
        double b2 = B[k + 2*k_size];
        double b3 = B[k + 3*k_size];

        c00 += a0 * b0;  c01 += a0 * b1;  c02 += a0 * b2;  c03 += a0 * b3;
        c10 += a1 * b0;  c11 += a1 * b1;  c12 += a1 * b2;  c13 += a1 * b3;
        c20 += a2 * b0;  c21 += a2 * b1;  c22 += a2 * b2;  c23 += a2 * b3;
        c30 += a3 * b0;  c31 += a3 * b1;  c32 += a3 * b2;  c33 += a3 * b3;
    }

    // Store results back
    C(0,0) = c00;  C(0,1) = c01;  C(0,2) = c02;  C(0,3) = c03;
    C(1,0) = c10;  C(1,1) = c11;  C(1,2) = c12;  C(1,3) = c13;
    C(2,0) = c20;  C(2,1) = c21;  C(2,2) = c22;  C(2,3) = c23;
    C(3,0) = c30;  C(3,1) = c31;  C(3,2) = c32;  C(3,3) = c33;
}

// Pack block of matrix A into contiguous memory
static inline void pack_a_k22(const double* __restrict__ A, double* __restrict__ A_packed, int LDA, int M, int K) {
    int mp = (M + MR - 1) & -MR;  // Round up to multiple of MR
    
    // Pack complete MR x K blocks
    for (int i = 0; i < M - MR + 1; i += MR) {
        for (int k = 0; k < K; k++) { // loop unrolled
            for (int ii = 0; ii < MR; ii++) { // was vectorized
                A_packed[ii + k*MR] = A(i + ii, k);
            }
        }
        A_packed += K * MR; // sinked after loop using last value computation
    }
    
    // Pack remaining rows if M is not multiple of MR
    if (M != mp) {
        for (int k = 0; k < K; k++) {
            int i;
            for (i = 0; i < M % MR; i++) { // was vectorized, length=2
                A_packed[i + k*MR] = A(M - M % MR + i, k);
            }
            for (; i < MR; i++) {
                A_packed[i + k*MR] = 0.0;
            }
        }
    }
}

// Pack block of matrix B into contiguous memory
static inline void pack_b_k22(const double* __restrict__ B, double* __restrict__ B_packed, int LDB, int K, int N) {
    int np = (N + NR - 1) & -NR;  // Round up to multiple of NR
    
    // Pack complete K x NR blocks
    for (int j = 0; j < N - NR + 1; j += NR) {
        for (int jj = 0; jj < NR; jj++) { 
#pragma omp simd
            for (int k = 0; k < K; k++) { 
                B_packed[k + jj*K] = B(k, j + jj);
            }
        }
        B_packed += K * NR;
    }
    
    // Pack remaining columns if N is not multiple of NR
    if (N != np) {
        int j;
        for (j = 0; j < N % NR; j++) { // vector dependence prevents vectorization
#pragma omp simd
            for (int k = 0; k < K; k++) {
                B_packed[k + j*K] = B(k, N - N % NR + j);
            }
        }
        for (; j < NR; j++) { // vectorization possible but seems inefficient
#pragma omp simd
            for (int k = 0; k < K; k++) {
                B_packed[k + j*K] = 0.0;
            }
        }
    }
}

// Process a block of matrices
static void macro_kernel_k22(int M, int N, int K,
                           double alpha,
                           const double *A_packed,
                           const double *B_packed,
                           double *C, int LDC) {
    for (int j = 0; j < N - NR + 1; j += NR) {
        for (int i = 0; i < M - MR + 1; i += MR) {
            micro_kernel_k22(K, alpha,
                           A_packed + i*K,
                           B_packed + j*K,
                           &C(i,j), LDC);
        }
        
        // Handle edge case for M
        if (M % MR != 0) {
            for (int k = 0; k < K; k++) {
                for (int i = M - M % MR; i < M; i++) {
                    double ai = alpha * A_packed[i + k*M];
                    for (int jj = 0; jj < NR; jj++) { // vector dependence prevents vectorization, unrolled by 4
                        C(i,j+jj) += ai * B_packed[k + jj*K];
                    }
                }
            }
        }
    }
    
    // Handle edge case for N
    if (N % NR != 0) {
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                double ai = alpha * A_packed[i + k*M];
                for (int j = N - N % NR; j < N; j++) { // vector dependence prevents vectorization, unrolled by 4
                    C(i,j) += ai * B_packed[k + j*K];
                }
            }
        }
    }
}

// Main DGEMM function
void mydgemm_cpu_v22(int M, int N, int K,
                     double alpha, const double *A, int LDA,
                     const double *B, int LDB,
                     double beta, double *C, int LDC) {
    
    // Handle special cases
    if (M == 0 || N == 0 || K == 0) return;
    if (alpha == 0.0) {
        scale_c_k22(C, M, N, LDC, beta);
        return;
    }

    // Scale matrix C by beta
    scale_c_k22(C, M, N, LDC, beta);
    
    // Calculate sizes for packed buffers
    int mp = (M + MR - 1) & -MR;  // Round up to multiple of MR
    int np = (N + NR - 1) & -NR;  // Round up to multiple of NR
    
    // Allocate packed buffers
    double *A_packed = (double *)aligned_alloc(64, mp * K * sizeof(double));
    double *B_packed = (double *)aligned_alloc(64, K * np * sizeof(double));
    
    if (!A_packed || !B_packed) {
        if (A_packed) free(A_packed);
        if (B_packed) free(B_packed);
        return;  // Memory allocation failed
    }
    
    // Block loops
    for (int n_start = 0; n_start < N; n_start += N_BLOCK) {
        int n_size = min(N_BLOCK, N - n_start);
        
        for (int k_start = 0; k_start < K; k_start += K_BLOCK) {
            int k_size = min(K_BLOCK, K - k_start);
            
            // Pack B block
            pack_b_k22(&B(k_start, n_start), B_packed, LDB, k_size, n_size);
            
            for (int m_start = 0; m_start < M; m_start += M_BLOCK) {
                int m_size = min(M_BLOCK, M - m_start);
                
                // Pack A block
                pack_a_k22(&A(m_start, k_start), A_packed, LDA, m_size, k_size);
                
                // Compute C block
                macro_kernel_k22(m_size, n_size, k_size, alpha,
                               A_packed, B_packed,
                               &C(m_start, n_start), LDC);
            }
        }
    }
    
    free(A_packed);
    free(B_packed);
}

#endif // KERNEL22_H
