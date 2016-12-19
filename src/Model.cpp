#include "Model.h"

#include <cstring>
#include "GEMM.h"

void computeA(Material const& material, double A[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES])
{
  memset(A, 0, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES*sizeof(double));
  A[0 * NUMBER_OF_QUANTITIES + 1] = material.K0;
  A[1 * NUMBER_OF_QUANTITIES + 0] = 1.0 / material.rho0;
}

void computeB(Material const& material, double B[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES])
{
  memset(B, 0, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES*sizeof(double));
  B[0 * NUMBER_OF_QUANTITIES + 2] = material.K0;
  B[2 * NUMBER_OF_QUANTITIES + 0] = 1.0 / material.rho0;
}

void rotateFluxSolver(  double        nx,
                        double        ny,
                        double const  fluxSolver[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                        double        rotatedFluxSolver[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] )
{
  double T[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] = {}; // zero initialisation
  double TT[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] = {}; // zero initialisation
  double tmp[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] = {}; // zero initialisation
  
  memset(rotatedFluxSolver, 0, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES*sizeof(double));
  
  T[0*NUMBER_OF_QUANTITIES + 0] = 1.0;
  T[1*NUMBER_OF_QUANTITIES + 1] = nx;
  T[1*NUMBER_OF_QUANTITIES + 2] = ny;
  T[2*NUMBER_OF_QUANTITIES + 1] = -ny;
  T[2*NUMBER_OF_QUANTITIES + 2] = nx;
  
  TT[0*NUMBER_OF_QUANTITIES + 0] = 1.0;
  TT[1*NUMBER_OF_QUANTITIES + 1] = nx;
  TT[1*NUMBER_OF_QUANTITIES + 2] = -ny;
  TT[2*NUMBER_OF_QUANTITIES + 1] = ny;
  TT[2*NUMBER_OF_QUANTITIES + 2] = nx;
  
  DGEMM(  NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
          1.0, T, NUMBER_OF_QUANTITIES,
          fluxSolver, NUMBER_OF_QUANTITIES,
          0.0, tmp, NUMBER_OF_QUANTITIES );
  
  DGEMM(  NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES,
          1.0, tmp, NUMBER_OF_QUANTITIES,
          TT, NUMBER_OF_QUANTITIES,
          0.0, rotatedFluxSolver, NUMBER_OF_QUANTITIES );
}

void computeAplus( Material const&  local,
                   Material const&  neighbour,
                   double           Aplus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] )
{
  memset(Aplus, 0, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES*sizeof(double));
  
  double cm = local.wavespeed();
  double cp = neighbour.wavespeed();
  double div1 = 1.0 / (local.K0 * cp + neighbour.K0 * cm);
  double div2 = div1 / local.rho0;
  Aplus[0*NUMBER_OF_QUANTITIES + 0] = local.K0 * cm * cp * div1;
  Aplus[0*NUMBER_OF_QUANTITIES + 1] = local.K0 * local.K0 * cp * div1;
  Aplus[1*NUMBER_OF_QUANTITIES + 0] = neighbour.K0 * cm * div2;
  Aplus[1*NUMBER_OF_QUANTITIES + 1] = local.K0 * neighbour.K0 * div2;
}

void computeAminus( Material const& local,
                    Material const& neighbour,
                    double          Aminus[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] )
{  
  memset(Aminus, 0, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES*sizeof(double));
  
  double cm = local.wavespeed();
  double cp = neighbour.wavespeed();
  double div1 = 1.0 / (local.K0 * cp + neighbour.K0 * cm);
  double div2 = div1 / local.rho0;

  Aminus[0*NUMBER_OF_QUANTITIES + 0] = -local.K0 * cm * cp * div1;
  Aminus[0*NUMBER_OF_QUANTITIES + 1] = local.K0 * neighbour.K0 * cm * div1;
  Aminus[1*NUMBER_OF_QUANTITIES + 0] = local.K0 * cp * div2;
  Aminus[1*NUMBER_OF_QUANTITIES + 1] = -local.K0 * neighbour.K0 * div2;
}
