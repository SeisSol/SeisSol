/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 * 
 **/

#include <PUML/PUML.h>
#include <PUML/Downward.h>
#include "LtsWeights.h"

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <Initializer/ParameterDB.h>

void seissol::initializers::time_stepping::LtsWeights::computeMaxTimesteps( PUML::TETPUML const&  mesh,
                                                                            double const*         lambda,
                                                                            double const*         mu,
                                                                            double const*         rho,
                                                                            double*               timestep ) {
  std::vector<PUML::TETPUML::cell_t> const& cells = mesh.cells();
	std::vector<PUML::TETPUML::vertex_t> const& vertices = mesh.vertices();

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    // Compute maximum wavespeed
    double pWaveVel = sqrt( (lambda[cell] + 2.0 * mu[cell]) / rho[cell] );
    
    // Compute insphere radius
    glm::dvec3 barycentre(0.,0.,0.);
    glm::dvec3 x[4];
    unsigned vertLids[4];
    PUML::Downward::vertices(mesh, cells[cell], vertLids);
    for (unsigned vtx = 0; vtx < 4; ++vtx) {
      for (unsigned d = 0; d < 3; ++d) {
        x[vtx][d] = vertices[ vertLids[vtx] ].coordinate()[d];
      }
    }

    double alpha = determinant(glm::dmat4(glm::dvec4(x[0], 1.0), glm::dvec4(x[1], 1.0), glm::dvec4(x[2], 1.0), glm::dvec4(x[3], 1.0)));
    double Nabc = length(cross(x[1]-x[0], x[2]-x[0]));
    double Nabd = length(cross(x[1]-x[0], x[3]-x[0]));
    double Nacd = length(cross(x[2]-x[0], x[3]-x[0]));
    double Nbcd = length(cross(x[2]-x[1], x[3]-x[1]));
    double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);
    
    // Compute maximum timestep (CFL=1)
    timestep[cell] = 2.0 * insphere / (pWaveVel * (2*CONVERGENCE_ORDER-1));
  }
}

int seissol::initializers::time_stepping::LtsWeights::getCluster( double    timestep,
                                                                  double    globalMinTimestep,
                                                                  unsigned  rate  ) {
  if (rate == 1) {
    return 0;
  }

  double upper;
  upper = rate * globalMinTimestep;

  int cluster = 0;
  while (upper <= timestep) {
    upper *= rate;
    ++cluster;
  }
  return cluster;
}

int seissol::initializers::time_stepping::LtsWeights::ipow(int x, int y) {
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while(--y) {
    result *= x;
  }
  return result;
}

void seissol::initializers::time_stepping::LtsWeights::computeWeights(PUML::TETPUML const& mesh) {
  std::vector<PUML::TETPUML::cell_t> const& cells = mesh.cells();
  int const* boundaryCond = mesh.cellData(1);

  double* rho = new double[cells.size()];
  double* mu = new double[cells.size()];
  double* lambda = new double[cells.size()];
  
  seissol::initializers::ElementBarycentreGeneratorPUML queryGen(mesh);  
  seissol::initializers::ParameterDB parameterDB;
  parameterDB.addParameter("rho", rho);
  parameterDB.addParameter("mu", mu);
  parameterDB.addParameter("lambda", lambda);
  parameterDB.evaluateModel(m_velocityModel, queryGen);
  
  double* timestep = new double[cells.size()];  
  computeMaxTimesteps(mesh, lambda, mu, rho, timestep);
    
  delete[] lambda;
  delete[] mu;
  delete[] rho;

  double localMinTimestep = *std::min_element(timestep, timestep + cells.size());
  double localMaxTimestep = *std::max_element(timestep, timestep + cells.size());
  double globalMinTimestep;
  double globalMaxTimestep;
#ifdef USE_MPI
  MPI_Allreduce(&localMinTimestep, &globalMinTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&localMaxTimestep, &globalMaxTimestep, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  globalMinTimestep = localMinTimestep;
  globalMaxTimestep = localMaxTimestep;
#endif
  
  delete[] m_vertexWeights;
  m_vertexWeights = new int[cells.size()];
  int maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, m_rate);
  int drToCellRatio = 1;
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int cluster = getCluster(timestep[cell], globalMinTimestep, m_rate);
    
    int dynamicRupture = 0;
    for (unsigned face = 0; face < 4; ++face) {
      dynamicRupture += ( (boundaryCond[cell] >> (face*8)) & 0xFF == 3) ? 1 : 0;
    }
    
    m_vertexWeights[cell] = (1 + drToCellRatio*dynamicRupture) * ipow(m_rate, maxCluster - cluster);
  }
  
  delete[] timestep;
}
