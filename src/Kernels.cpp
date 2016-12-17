#include "Kernels.h"
#include "GEMM.h"
#include "constants.h"
#include "GlobalMatrices.h"
#include "Model.h"

void computeAder( double                  timestep,
                  GlobalConstants const&  globals,
                  Material const&         material,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated )
{
  double tmp[NUMBER_OF_DOFS] = {}; // zero initialisation
  double derivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] = {}; // zero initialisation
  double A[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  double B[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  
  double factor = timestep;

  for (unsigned dof = 0; dof < NUMBER_OF_DOFS; ++dof) {
    derivatives[0][dof] = degreesOfFreedom[dof];
    timeIntegrated[dof] = factor * degreesOfFreedom[dof];
  }
  
  computeA(material, A);
  computeB(material, B);
  
  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {  
    // tmp = Kxi^T * degreesOfFreedom
    DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
            1.0, GlobalMatrices::KxiT, NUMBER_OF_BASIS_FUNCTIONS,
            derivatives[der-1], NUMBER_OF_BASIS_FUNCTIONS,
            0.0, tmp, NUMBER_OF_BASIS_FUNCTIONS );
    
    // derivatives[der] = -1/hx * tmp * A
    DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
            -1.0 / globals.hx, tmp, NUMBER_OF_BASIS_FUNCTIONS,
            A, NUMBER_OF_QUANTITIES,
            1.0, derivatives[der], NUMBER_OF_BASIS_FUNCTIONS );
    
    // tmp = Keta^T * degreesOfFreedom
    DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
            1.0, GlobalMatrices::KetaT, NUMBER_OF_BASIS_FUNCTIONS,
            derivatives[der-1], NUMBER_OF_BASIS_FUNCTIONS,
            0.0, tmp, NUMBER_OF_BASIS_FUNCTIONS );
    
    // derivatives[der] += -1/hy * tmp * B
    DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
            -1.0 / globals.hy, tmp, NUMBER_OF_BASIS_FUNCTIONS,
            B, NUMBER_OF_QUANTITIES,
            1.0, derivatives[der], NUMBER_OF_BASIS_FUNCTIONS );

    factor *= timestep / (der + 1);
    for (unsigned dof = 0; dof < NUMBER_OF_DOFS; ++dof) {
      timeIntegrated[dof] += factor * derivatives[der][dof];
    }
  }
}

void computeVolumeIntegral( GlobalConstants const&  globals,
                            Material const&         material,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom )
{
  double A[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  double B[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  double tmp[NUMBER_OF_DOFS] = {}; // zero initialisation
  
  computeA(material, A);
  computeB(material, B);
  
  // Computes tmp = Kxi * timeIntegrated
  DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
          1.0, GlobalMatrices::Kxi, NUMBER_OF_BASIS_FUNCTIONS,
          timeIntegrated, NUMBER_OF_BASIS_FUNCTIONS,
          0.0, tmp, NUMBER_OF_BASIS_FUNCTIONS );
  
  // Computes degreesOfFreedom += 1/hx tmp * A
  DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
          1.0 / globals.hx, tmp, NUMBER_OF_BASIS_FUNCTIONS,
          A, NUMBER_OF_QUANTITIES,
          1.0, degreesOfFreedom, NUMBER_OF_BASIS_FUNCTIONS );
  
  // Computes tmp = Keta * timeIntegrated
  DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
          1.0, GlobalMatrices::Keta, NUMBER_OF_BASIS_FUNCTIONS,
          timeIntegrated, NUMBER_OF_BASIS_FUNCTIONS,
          0.0, tmp, NUMBER_OF_BASIS_FUNCTIONS );
  
  // Computes degreesOfFreedom += 1/hy tmp * B
  DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
          1.0 / globals.hy, tmp, NUMBER_OF_BASIS_FUNCTIONS,
          B, NUMBER_OF_QUANTITIES,
          1.0, degreesOfFreedom, NUMBER_OF_BASIS_FUNCTIONS );
}

void computeFlux( double                  factor,
                  double const            fluxMatrix[NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_BASIS_FUNCTIONS],
                  double const            rotatedFluxSolver[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                  DegreesOfFreedom const& timeIntegrated,
                  DegreesOfFreedom        degreesOfFreedom )
{
  double tmp[NUMBER_OF_DOFS] = {}; // zero initialisation
  
  // Computes tmp = fluxMatrix * timeIntegrated
  DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_BASIS_FUNCTIONS,
          1.0, fluxMatrix, NUMBER_OF_BASIS_FUNCTIONS,
          timeIntegrated, NUMBER_OF_BASIS_FUNCTIONS,
          0.0, tmp, NUMBER_OF_BASIS_FUNCTIONS );
  
  // Computes degreesOfFreedom += factor * tmp * rotatedFluxSolver
  DGEMM(  NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
          factor, tmp, NUMBER_OF_BASIS_FUNCTIONS,
          rotatedFluxSolver, NUMBER_OF_QUANTITIES,
          1.0, degreesOfFreedom, NUMBER_OF_BASIS_FUNCTIONS );
}
