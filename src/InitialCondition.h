#ifndef INITIALCONDITION_H_
#define INITIALCONDITION_H_

#include "typedefs.h"
#include "Grid.h"

void initialCondition(  GlobalConstants const& globals,
                        Grid<Material>& materialGrid,
                        Grid<DegreesOfFreedom>& degreesOfFreedomGrid  );

void L2error( double time,
              GlobalConstants const& globals,
              Grid<Material>& materialGrid,
              Grid<DegreesOfFreedom>& degreesOfFreedomGrid,
              double l2error[NUMBER_OF_QUANTITIES]  );

void initSourcetermPhi(double xi, double eta, SourceTerm& sourceterm);

#endif // INITIALCONDITION_H_
