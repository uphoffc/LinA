#ifndef KERNELS_H_
#define KERNELS_H_

#include "typedefs.h"

void computeAder( double                  timestep,
                  GlobalConstants const&  globals,
                  Material const&         material,
                  DegreesOfFreedom const& degreesOfFreedom,
                  DegreesOfFreedom&       timeIntegrated );
                  
void computeVolumeIntegral( GlobalConstants const&  globals,
                            Material const&         material,
                            DegreesOfFreedom const& timeIntegrated,
                            DegreesOfFreedom&       degreesOfFreedom );

void computeFlux( double                  factor,
                  double const            fluxMatrix[NUMBER_OF_BASIS_FUNCTIONS*NUMBER_OF_BASIS_FUNCTIONS],
                  double const            rotatedFluxSolver[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES],
                  DegreesOfFreedom const& timeIntegrated,
                  DegreesOfFreedom        degreesOfFreedom );

#endif // KERNELS_H_
