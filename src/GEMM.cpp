#include "GEMM.h"

void DGEMM(unsigned M, unsigned N, unsigned K, double alpha, double const* A, unsigned ldA, double const* B, unsigned ldB, double beta, double* C, unsigned ldC) {
  for (unsigned j = 0; j < N; ++j) {
    for (unsigned i = 0; i < M; ++i) {
      double cij = 0.0;
      for (unsigned k = 0; k < K; ++k) {
        cij += A[k*ldA + i] * B[j*ldB + k];
      }
      C[j*ldC + i] = alpha * cij + beta * C[j*ldC + i];
    }
  }
}
