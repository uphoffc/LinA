#ifndef WAVEFIELDWRITER_H_
#define WAVEFIELDWRITER_H_

#include <string>
#include <fstream>
#include "typedefs.h"
#include "Grid.h"

class WaveFieldWriter {
public:
  WaveFieldWriter(std::string const& baseName, GlobalConstants const& globals, double interval, int pointsPerDim);
  ~WaveFieldWriter();
  
  void writeTimestep(double time, Grid<DegreesOfFreedom>& degreesOfFreedomGrid, bool forceWrite = false);
private:
  unsigned      m_step;
  std::string   m_dirName;
  std::string   m_baseName;
  std::ofstream m_xdmf;
  float*        m_pressure;
  float*        m_uvel;
  float*        m_vvel;
  double        m_interval;
  double        m_lastTime;
  double*       m_subsampleMatrix;
  double*       m_subsamples;
  int           m_pointsPerDim;
};


#endif // WAVEFIELDWRITER_H_
