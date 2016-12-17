#ifndef GEMM_H_
#define GEMM_H_

/// Generalized matrix-matrix multiplication for column-major storage
void DGEMM(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC);

#endif // GEMM_H_
