#include "WaveFieldWriter.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <limits>

WaveFieldWriter::WaveFieldWriter(std::string const& baseName, GlobalConstants const& globals, double interval)
  : m_step(0), m_interval(interval), m_lastTime(-std::numeric_limits<double>::max())
{
  if (!baseName.empty()) {
    m_xdmf.open((baseName + ".xdmf").c_str());
    m_xdmf  << "<?xml version=\"1.0\" ?>" << std::endl
            << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">" << std::endl
            << "<Xdmf Version=\"2.0\">" << std::endl
            << "  <Domain>" << std::endl
            << "    <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << globals.Y+1 << " " << globals.X+1 << "\"/>" << std::endl
            << "    <Geometry GeometryType=\"ORIGIN_DXDY\">" << std::endl
            << "      <DataItem Format=\"XML\" Dimensions=\"2\">0.0 0.0</DataItem>" << std::endl
            << "      <DataItem Format=\"XML\" Dimensions=\"2\">" << globals.hy << " " << globals.hx << "</DataItem>" << std::endl
            << "    </Geometry>" << std::endl
            << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
            
    m_pressure = new float[globals.X*globals.Y];
    m_uvel = new float[globals.X*globals.Y];
    m_vvel = new float[globals.X*globals.Y];
    
    std::size_t lastFound = 0;
    std::size_t found;
    while ((found = baseName.find("/", lastFound+1)) != std::string::npos) {
      lastFound = found;
    }
    if (lastFound > 0) {
      ++lastFound;
    }
    m_dirName = baseName.substr(0, lastFound);
    m_baseName = baseName.substr(lastFound);
  }
}

WaveFieldWriter::~WaveFieldWriter()
{
  if (!m_baseName.empty()) {
    m_xdmf  << "    </Grid>" << std::endl
            << "  </Domain>" << std::endl
            << "</Xdmf>" << std::endl;
    m_xdmf.close();

    delete[] m_pressure;
    delete[] m_uvel;
    delete[] m_vvel;
  }
}

void WaveFieldWriter::writeTimestep(double time, Grid<DegreesOfFreedom>& degreesOfFreedomGrid, bool forceWrite)
{
  if (!m_baseName.empty() && (time >= m_lastTime + m_interval || forceWrite)) {
    m_lastTime = time;
    
    std::stringstream pressureFileName, uvelFileName, vvelFileName;
    pressureFileName << m_baseName << "_pressure" << m_step << ".bin";
    uvelFileName << m_baseName << "_u" << m_step << ".bin";
    vvelFileName << m_baseName << "_v" << m_step << ".bin";
    
    m_xdmf  << "      <Grid Name=\"step_" << m_step << "\" GridType=\"Uniform\">" << std::setw(0) << std::endl
            << "        <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>" << std::endl
            << "        <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>" << std::endl
            << "        <Time Value=\"" << time << "\"/>" << std::endl
            << "        <Attribute Name=\"pressure\" Center=\"Cell\">" << std::endl
            << "          <DataItem Format=\"Binary\" DataType=\"Float\" Precision=\"4\" Dimensions=\"" << degreesOfFreedomGrid.Y() << " " << degreesOfFreedomGrid.X() << "\">" << std::endl
            << "            " << pressureFileName.str() << std::endl
            << "          </DataItem>" << std::endl
            << "        </Attribute>" << std::endl
            << "        <Attribute Name=\"u\" Center=\"Cell\">" << std::endl
            << "          <DataItem Format=\"Binary\" DataType=\"Float\" Precision=\"4\" Dimensions=\"" << degreesOfFreedomGrid.Y() << " " << degreesOfFreedomGrid.X() << "\">" << std::endl
            << "            " << uvelFileName.str() << std::endl
            << "          </DataItem>" << std::endl
            << "       </Attribute>" << std::endl
            << "        <Attribute Name=\"v\" Center=\"Cell\">" << std::endl
            << "          <DataItem Format=\"Binary\" DataType=\"Float\" Precision=\"4\" Dimensions=\"" << degreesOfFreedomGrid.Y() << " " << degreesOfFreedomGrid.X() << "\">" << std::endl
            << "            " << vvelFileName.str() << std::endl
            << "          </DataItem>" << std::endl
            << "        </Attribute>" << std::endl
            << "      </Grid>" << std::endl;

    for (int y = 0; y < degreesOfFreedomGrid.Y(); ++y) {
      for (int x = 0; x < degreesOfFreedomGrid.X(); ++x) {
        // 0-th degree of freedom is equal to cell-average due to choice of basis functions
        DegreesOfFreedom& degreesOfFreedom = degreesOfFreedomGrid.get(x, y);
        m_pressure[y*degreesOfFreedomGrid.X() + x] = degreesOfFreedom[0 * NUMBER_OF_BASIS_FUNCTIONS + 0];
        m_uvel[y*degreesOfFreedomGrid.X() + x] = degreesOfFreedom[1 * NUMBER_OF_BASIS_FUNCTIONS + 0];
        m_vvel[y*degreesOfFreedomGrid.X() + x] = degreesOfFreedom[2 * NUMBER_OF_BASIS_FUNCTIONS + 0];
      }
    }
    
    FILE* pressureFile = fopen((m_dirName + pressureFileName.str()).c_str(), "wb");
    fwrite(m_pressure, sizeof(float), degreesOfFreedomGrid.X()*degreesOfFreedomGrid.Y(), pressureFile);
    fclose(pressureFile);
    
    FILE* uFile = fopen((m_dirName + uvelFileName.str()).c_str(), "wb");
    fwrite(m_uvel, sizeof(float), degreesOfFreedomGrid.X()*degreesOfFreedomGrid.Y(), uFile);
    fclose(uFile);
    
    FILE* vFile = fopen((m_dirName + vvelFileName.str()).c_str(), "wb");
    fwrite(m_vvel, sizeof(float), degreesOfFreedomGrid.X()*degreesOfFreedomGrid.Y(), vFile);
    fclose(vFile);
    
    ++m_step;
  }
}
