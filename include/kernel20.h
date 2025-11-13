#ifndef KERNEL20_H
#define KERNEL20_H

#include <stdlib.h>

// Access macros for matrices in column-major format
#define A(i,j) A[(i)+(j)*LDA]
#define B(i,j) B[(i)+(j)*LDB]
#define C(i,j) C[(i)+(j)*LDC]

// Cache blocking parameters - adjust based on your CPU's cache sizes
#define M_BLOCKING 192  // Block size for matrix rows
#define N_BLOCKING 2048 // Block size for matrix columns
#define K_BLOCKING 384  // Block size for inner dimension

// Helper function to get minimum of two values
static inline int min(int a, int b) {
    return (a < b) ? a : b;
}

// Scale matrix C by beta
void scale_matrix_c_k20(double *C, int M, int N, int LDC, double beta) {
    if (beta == 0.0) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) {
                C(i,j) = 0.0;
            }
        }
    } else if (beta != 1.0) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < M; i++) {
                C(i,j) *= beta;
            }
        }
    }
}

// Pack matrix A into contiguous memory
// src is in column-major format
// dst will be packed for sequential access in block multiplication
void pack_matrix_a_k20(double *src, double *dst, int LDA, int rows, int cols) {
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            dst[i + j*rows] = src[i + j*LDA];
        }
    }
}

// Pack matrix B into contiguous memory
// src is in column-major format
// dst will be packed for sequential access in block multiplication
void pack_matrix_b_k20(double *src, double *dst, int LDB, int rows, int cols) {
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            dst[i*cols + j] = src[i + j*LDB];
        }
    }
}

// Compute C += alpha * A * B for a single block
// A and B are packed, C is in column-major format
void compute_block_k20(const double *A, const double *B, double *C, int LDC,
                      int M, int N, int K, double alpha) {
    // Process the matrix C in 4x4 blocks for better register reuse
    for (int j = 0; j < N; j += 4) {
        int N_block = min(4, N - j);
        
        for (int i = 0; i < M; i += 4) {
            int M_block = min(4, M - i);
            
            // Initialize accumulators for 4x4 block of C
            double c[4][4] = {{0.0}};
            
            // Compute contribution from A and B
            for (int k = 0; k < K; k++) {
                // Load elements from packed A
                double a[4];
                for (int ii = 0; ii < M_block; ii++) {
                    a[ii] = alpha * A[i + ii + k*M];
                }
                
                // Load elements from packed B
                double b[4];
                for (int jj = 0; jj < N_block; jj++) {
                    b[jj] = B[k*N + j + jj];
                }
                
                // Outer product update of the C block
                for (int ii = 0; ii < M_block; ii++) {
                    for (int jj = 0; jj < N_block; jj++) {
                        c[ii][jj] += a[ii] * b[jj];
                    }
                }
            }
            
            // Store results back to C
            for (int ii = 0; ii < M_block; ii++) {
                for (int jj = 0; jj < N_block; jj++) {
                    C(i + ii, j + jj) += c[ii][jj];
                }
            }
        }
    }
}

// Main DGEMM function
// Computes C = beta * C + alpha * A * B
// Matrices are in column-major format
void mydgemm_cpu_v20(int M, int N, int K,
                     double alpha, const double *A, int LDA,
                     const double *B, int LDB,
                     double beta, double *C, int LDC) {
    
    // Scale matrix C by beta
    scale_matrix_c_k20(C, M, N, LDC, beta);
    
    // Allocate packed buffers with alignment for better memory access
    double *b_buffer = (double *)aligned_alloc(64, K_BLOCKING * N_BLOCKING * sizeof(double));
    double *a_buffer = (double *)aligned_alloc(64, K_BLOCKING * M_BLOCKING * sizeof(double));
    
    // Block loops
    for (int n_start = 0; n_start < N; n_start += N_BLOCKING) {
        int n_size = min(N_BLOCKING, N - n_start);
        
        for (int k_start = 0; k_start < K; k_start += K_BLOCKING) {
            int k_size = min(K_BLOCKING, K - k_start);
            
            // Pack B block
            pack_matrix_b_k20(B + k_start + n_start*LDB,
                             b_buffer, LDB,
                             k_size, n_size);
            
            for (int m_start = 0; m_start < M; m_start += M_BLOCKING) {
                int m_size = min(M_BLOCKING, M - m_start);
                
                // Pack A block
                pack_matrix_a_k20(A + m_start + k_start*LDA,
                                 a_buffer, LDA,
                                 m_size, k_size);
                
                // Compute C block
                compute_block_k20(a_buffer, b_buffer,
                                C + m_start + n_start*LDC, LDC,
                                m_size, n_size, k_size, alpha);
            }
        }
    }
    
    // Free packed buffers
    free(a_buffer);
    free(b_buffer);
}

#endif // KERNEL20_H