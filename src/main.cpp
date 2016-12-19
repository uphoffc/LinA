#include <cmath>
#include <iostream>

#include "tclap/CmdLine.h"
#include "typedefs.h"
#include "Simulator.h"
#include "Grid.h"
#include "InitialCondition.h"

void initScenario0(GlobalConstants& globals, Grid<Material>& materialGrid, Grid<DegreesOfFreedom>& degreesOfFreedomGrid)
{    
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      Material& material = materialGrid.get(x, y);
      material.rho0 = 1.;
      material.K0 = 4.;
    }
  }

  initialCondition(globals, materialGrid, degreesOfFreedomGrid);
}

void initScenario1(GlobalConstants& globals, Grid<Material>& materialGrid, Grid<DegreesOfFreedom>& degreesOfFreedomGrid)
{
  double checkerWidth = 0.25;

  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      Material& material = materialGrid.get(x, y);
      int matId = static_cast<int>(x*globals.hx/checkerWidth) % 2 ^ static_cast<int>(y*globals.hy/checkerWidth) % 2;
      if (matId == 0) {
        material.rho0 = 1.;
        material.K0 = 2.;
      } else {
        material.rho0 = 2.;
        material.K0 = 0.5;
      }
    }
  }

  initialCondition(globals, materialGrid, degreesOfFreedomGrid);
}

double sourceFunctionAntiderivative(double time)
{
  return sin(time);
}

void initSourceTerm23(GlobalConstants& globals, SourceTerm& sourceterm)
{
  sourceterm.quantity = 0; // pressure source
  double xs = 0.5;
  double ys = 0.5;
  sourceterm.x = static_cast<int>(xs / (globals.hx));
  sourceterm.y = static_cast<int>(ys / (globals.hy));
  double xi = (xs - sourceterm.x*globals.hx) / globals.hx;
  double eta = (ys - sourceterm.y*globals.hy) / globals.hy;
  
  initSourcetermPhi(xi, eta, sourceterm);
  
  sourceterm.antiderivative = sourceFunctionAntiderivative;  
}

void initScenario2(GlobalConstants& globals, Grid<Material>& materialGrid, Grid<DegreesOfFreedom>& degreesOfFreedomGrid, SourceTerm& sourceterm)
{
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      Material& material = materialGrid.get(x, y);
      material.rho0 = 1.;
      material.K0 = 2.;
    }
  }
  
  initSourceTerm23(globals, sourceterm);
}

void initScenario3(GlobalConstants& globals, Grid<Material>& materialGrid, Grid<DegreesOfFreedom>& degreesOfFreedomGrid, SourceTerm& sourceterm)
{
  for (int y = 0; y < globals.Y; ++y) {
    for (int x = 0; x < globals.X; ++x) {
      Material& material = materialGrid.get(x, y);
      int matId;
      double xp = x*globals.hx;
      double yp = y*globals.hy;
      matId = (xp >= 0.25 && xp <= 0.75 && yp >= 0.25 && yp <= 0.75) ? 0 : 1;
      if (matId == 0) {
        material.rho0 = 1.;
        material.K0 = 2.;
      } else {
        material.rho0 = 2.;
        material.K0 = 0.5;
      }
    }
  }
  
  initSourceTerm23(globals, sourceterm);
}

int main(int argc, char** argv)
{
  int scenario;
  double wfwInterval;
  std::string wfwBasename;
  GlobalConstants globals;
  
  try {
    TCLAP::CmdLine cmd("ADER-DG for linear acoustics.", ' ', "0.1");
    TCLAP::ValueArg<int> scenarioArg("s", "scenario", "Scenario. 0=Convergence test. 1=Checkerboard.", true, 0, "int");
    TCLAP::ValueArg<int> XArg("x", "x-number-of-cells", "Number of cells in x direction.", true, 0, "int");
    TCLAP::ValueArg<int> YArg("y", "y-number-of-cells", "Number of cells in y direction.", true, 0, "int");
    TCLAP::ValueArg<std::string> basenameArg("o", "output", "Basename of wavefield writer output. Leave empty for no output.", false, "", "string");
    TCLAP::ValueArg<double> intervalArg("i", "interval", "Time interval of wavefield writer.", false, 0.1, "double");
    TCLAP::ValueArg<double> timeArg("t", "end-time", "Final simulation time.", false, 0.5, "double");
    cmd.add(scenarioArg);
    cmd.add(XArg);
    cmd.add(YArg);
    cmd.add(basenameArg);
    cmd.add(intervalArg);
    cmd.add(timeArg);
    
    cmd.parse(argc, argv);
    
    scenario = scenarioArg.getValue();
    globals.X = XArg.getValue();
    globals.Y = YArg.getValue();
    wfwBasename = basenameArg.getValue();
    wfwInterval = intervalArg.getValue();
    globals.endTime = timeArg.getValue();
    
    if (scenario < 0 || scenario > 3) {
      std::cerr << "Unknown scenario." << std::endl;
      return -1;
    }
    if (globals.X < 0 || globals.Y < 0) {
      std::cerr << "X or Y smaller than 0. Does not make sense." << std::endl;
      return -1;
    }
  } catch (TCLAP::ArgException &e) {
    std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
    return -1;
  }
  
  globals.hx = 1. / globals.X;
  globals.hy = 1. / globals.Y;
  
  Grid<DegreesOfFreedom> degreesOfFreedomGrid(globals.X, globals.Y);
  Grid<Material> materialGrid(globals.X, globals.Y);
  SourceTerm sourceterm;
  
  switch (scenario) {
    case 0:
      initScenario0(globals, materialGrid, degreesOfFreedomGrid);
      break;
    case 1:
      initScenario1(globals, materialGrid, degreesOfFreedomGrid);
      break;
    case 2:
      initScenario2(globals, materialGrid, degreesOfFreedomGrid, sourceterm);
      break;
    case 3:
      initScenario3(globals, materialGrid, degreesOfFreedomGrid, sourceterm);
      break;
    default:
      break;
  }
  
  globals.maxTimestep = determineTimestep(globals.hx, globals.hy, materialGrid);
  
  WaveFieldWriter waveFieldWriter(wfwBasename, globals, wfwInterval, static_cast<int>(ceil( sqrt(NUMBER_OF_BASIS_FUNCTIONS) )));

  int steps = simulate(globals, materialGrid, degreesOfFreedomGrid, waveFieldWriter, sourceterm);
  
  if (scenario == 0) {
    double l2error[NUMBER_OF_QUANTITIES];
    L2error(globals.endTime, globals, materialGrid, degreesOfFreedomGrid, l2error);
    std::cout << "L2 error analysis" << std::endl << "=================" << std::endl;
    std::cout << "Pressue (p):    " << l2error[0] << std::endl;
    std::cout << "X-Velocity (u): " << l2error[1] << std::endl;
    std::cout << "Y-Velocity (v): " << l2error[2] << std::endl;
  }
  
  std::cout << "Total number of timesteps: " << steps << std::endl;
  
  return 0;
}
