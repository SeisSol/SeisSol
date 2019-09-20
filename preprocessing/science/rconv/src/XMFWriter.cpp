/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "XMFWriter.h"
#include <fstream>
#include <math.h>


void writeSlip(std::vector<SRFPointSource> const& sources, std::ofstream& xdmfFile)
{  
  xdmfFile << "      <Attribute Name=\"Slip path length (m)\" Center=\"Node\">" << std::endl
           << "        <DataItem Dimensions=\"" << sources.size() << "\" DataType=\"Float\" Precision=\"8\" Format=\"XML\">" << std::endl;
  for (std::vector<SRFPointSource>::const_iterator source = sources.begin(); source != sources.end(); ++source) {
    std::size_t maxSteps = std::max(std::max(source->slipRate[0].size(), source->slipRate[1].size()), source->slipRate[2].size());
    double slip = 0.0;
    double lastAbsSlipRate = 0.0;
    // Slip path integration with trapezoidal rule
    for (std::size_t step = 0; step < maxSteps; ++step) {
      double sr[3];
      for (int d = 0; d < 3; ++d) {
        sr[d] = (step < source->slipRate[d].size()) ? source->slipRate[d][step] : 0.0;
      }
      double absSlipRate = sqrt(sr[0]*sr[0] + sr[1]*sr[1] + sr[2]*sr[2]);
      slip += source->dt * (absSlipRate + lastAbsSlipRate) / 2.0;
      lastAbsSlipRate = absSlipRate;
    }
    // Convert to meter and write
    xdmfFile << slip / 100.0 << std::endl;
  }
  xdmfFile << "        </DataItem>" << std::endl
           << "      </Attribute>" << std::endl;
}

void writeXMF(char const* filename, std::vector<SRFPointSource> const& sources, Map const& map)
{  
  std::ofstream xdmfFile;
  xdmfFile.open(filename);
  
  xdmfFile << "<?xml version=\"1.0\" ?>" << std::endl
           << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
           << "<Xdmf>" << std::endl
           << "  <Domain>" << std::endl
           << "    <Topology Type=\"Polyvertex\" NumberOfElements=\"" << sources.size() << "\"/>" << std::endl
           << "    <Geometry Type=\"XYZ\">" << std::endl
           << "      <DataItem Dimensions=\"" << sources.size() << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">" << std::endl;

  for (std::vector<SRFPointSource>::const_iterator source = sources.begin(); source != sources.end(); ++source) {
    double x, y, z;
#ifdef noproj
    x = source->longitude;
    y = source->latitude;
    z = source->depth;
#else
    map.map(source->longitude, source->latitude, source->depth, &x, &y, &z);
#endif
    xdmfFile << x << " " << y << " " << z << std::endl;
  }

  xdmfFile << "      </DataItem>" << std::endl
           << "    </Geometry>" << std::endl
           << "    <Grid Name=\"Fault\" GridType=\"Uniform\">" << std::endl
           << "      <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>" << std::endl
           << "      <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>" << std::endl;

  writeSlip(sources, xdmfFile);

  xdmfFile << "    </Grid>" << std::endl
           << "  </Domain>" << std::endl
           << "</Xdmf>" << std::endl;
           
  xdmfFile.close();
}
