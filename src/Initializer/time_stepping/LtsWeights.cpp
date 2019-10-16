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
#include <Eigen/Eigenvalues>
#include <Kernels/precision.hpp>

#include <PUML/PUML.h>
#include <PUML/Downward.h>
#include <PUML/Upward.h>
#include "LtsWeights.h"

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <Initializer/ParameterDB.h>
#include <Parallel/MPI.h>

#include <generated_code/tensor.h>
#include <generated_code/init.h>

class FaceSorter {
private:
	std::vector<PUML::TETPUML::face_t> const& m_faces;

public:
	FaceSorter(std::vector<PUML::TETPUML::face_t> const& faces) : m_faces(faces) {}

	bool operator()(unsigned int a, unsigned int b) const {
		return m_faces[a].gid() < m_faces[b].gid();
	}
};

void seissol::initializers::time_stepping::LtsWeights::computeMaxTimesteps( PUML::TETPUML const&  mesh,
                                                                            double const*         pWaveVel,
                                                                            double*               timestep ) {
  std::vector<PUML::TETPUML::cell_t> const& cells = mesh.cells();
  std::vector<PUML::TETPUML::vertex_t> const& vertices = mesh.vertices();

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
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
    timestep[cell] = 2.0 * insphere / (pWaveVel[cell] * (2*CONVERGENCE_ORDER-1));
  }
}

int seissol::initializers::time_stepping::LtsWeights::getCluster( double    timestep,
                                                                  double    globalMinTimestep,
                                                                  unsigned  rate ) {
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

int seissol::initializers::time_stepping::LtsWeights::getBoundaryCondition( int const* boundaryCond,
                                                                            unsigned cell,
                                                                            unsigned face ) {
  int bcCurrentFace = ((boundaryCond[cell] >> (face*8)) & 0xFF);
  if (bcCurrentFace > 6) {
     bcCurrentFace = 3;
  }
  return bcCurrentFace;
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
  logInfo(seissol::MPI::mpi.rank()) << "Computing LTS weights.";

  std::vector<PUML::TETPUML::cell_t> const& cells = mesh.cells();
  int const* boundaryCond = mesh.cellData(1);

  double* pWaveVel = new double[cells.size()];
  
  seissol::initializers::ElementBarycentreGeneratorPUML queryGen(mesh);  
  seissol::initializers::ParameterDB parameterDB;
#ifndef USE_ANISOTROPIC
  double* rho = new double[cells.size()];
  double* mu = new double[cells.size()];
  double* lambda = new double[cells.size()];
  parameterDB.addParameter("rho", rho);
  parameterDB.addParameter("mu", mu);
  parameterDB.addParameter("lambda", lambda);
  parameterDB.evaluateModel(m_velocityModel, queryGen);
  for(unsigned cell = 0; cell < cells.size(); ++cell) {
    pWaveVel[cell] = sqrt( (lambda[cell] + 2.0 * mu[cell]) / rho[cell] );
  }
  delete[] rho;
  delete[] mu;
  delete[] lambda;
#else
  double* rho = new double[cells.size()];
  double* c11 = new double[cells.size()];
  double* c12 = new double[cells.size()];
  double* c13 = new double[cells.size()];
  double* c14 = new double[cells.size()];
  double* c15 = new double[cells.size()];
  double* c16 = new double[cells.size()];
  double* c22 = new double[cells.size()];
  double* c23 = new double[cells.size()];
  double* c24 = new double[cells.size()];
  double* c25 = new double[cells.size()];
  double* c26 = new double[cells.size()];
  double* c33 = new double[cells.size()];
  double* c34 = new double[cells.size()];
  double* c35 = new double[cells.size()];
  double* c36 = new double[cells.size()];
  double* c44 = new double[cells.size()];
  double* c45 = new double[cells.size()];
  double* c46 = new double[cells.size()];
  double* c55 = new double[cells.size()];
  double* c56 = new double[cells.size()];
  double* c66 = new double[cells.size()];
  parameterDB.addParameter("rho", rho);
  parameterDB.addParameter("c11", c11);
  parameterDB.addParameter("c12", c12);
  parameterDB.addParameter("c13", c13);
  parameterDB.addParameter("c14", c14);
  parameterDB.addParameter("c15", c15);
  parameterDB.addParameter("c16", c16);
  parameterDB.addParameter("c22", c22);
  parameterDB.addParameter("c23", c23);
  parameterDB.addParameter("c24", c24);
  parameterDB.addParameter("c25", c25);
  parameterDB.addParameter("c26", c26);
  parameterDB.addParameter("c33", c33);
  parameterDB.addParameter("c34", c34);
  parameterDB.addParameter("c35", c35);
  parameterDB.addParameter("c36", c36);
  parameterDB.addParameter("c44", c44);
  parameterDB.addParameter("c45", c45);
  parameterDB.addParameter("c46", c46);
  parameterDB.addParameter("c55", c55);
  parameterDB.addParameter("c56", c56);
  parameterDB.addParameter("c66", c66);
  parameterDB.evaluateModel(m_velocityModel, queryGen);
 
  using Matrix33 = Eigen::Matrix<real, 3, 3, Eigen::ColMajor>;
  Eigen::SelfAdjointEigenSolver<Matrix33> saes; 

  double samplingDirectionsData[seissol::tensor::samplingDirections::Size];
  std::copy_n(init::samplingDirections::Values,
              seissol::tensor::samplingDirections::Size,
              samplingDirectionsData);
  auto samplingDirections = init::samplingDirections::view::create(samplingDirectionsData);

  //compute elementwise wavespeeds
  for(unsigned cell = 0; cell < cells.size(); ++cell)
  { 
    real maxEv = 0;
    for(unsigned j = 0; j < 200; ++j)
    {
      real n0 = samplingDirections(j, 0);
      real n1 = samplingDirections(j, 1);
      real n2 = samplingDirections(j, 2);
      real M11 = n0 * c11[cell] * n0 + n0 * c16[cell] * n1 + n0 * c15[cell] * n2 + n1 * c16[cell] * n0 + n1 * c66[cell] * n1 + n1 * c56[cell] * n2 + n2 * c15[cell] * n0 + n2 * c56[cell] * n1 + n2 * c55[cell] * n2;
      real M22 = n0 * c66[cell] * n0 + n0 * c26[cell] * n1 + n0 * c46[cell] * n2 + n1 * c26[cell] * n0 + n1 * c22[cell] * n1 + n1 * c24[cell] * n2 + n2 * c46[cell] * n0 + n2 * c24[cell] * n1 + n2 * c44[cell] * n2; 
      real M33 = n0 * c55[cell] * n0 + n0 * c45[cell] * n1 + n0 * c35[cell] * n2 + n1 * c45[cell] * n0 + n1 * c44[cell] * n1 + n1 * c34[cell] * n2 + n2 * c35[cell] * n0 + n2 * c34[cell] * n1 + n2 * c33[cell] * n2;
      real M12 = n0 * c16[cell] * n0 + n0 * c12[cell] * n1 + n0 * c14[cell] * n2 + n1 * c66[cell] * n0 + n1 * c26[cell] * n1 + n1 * c46[cell] * n2 + n2 * c56[cell] * n0 + n2 * c25[cell] * n1 + n2 * c45[cell] * n2; 
      real M23 = n0 * c56[cell] * n0 + n0 * c46[cell] * n1 + n0 * c36[cell] * n2 + n1 * c25[cell] * n0 + n1 * c24[cell] * n1 + n1 * c23[cell] * n2 + n2 * c45[cell] * n0 + n2 * c44[cell] * n1 + n2 * c34[cell] * n2; 
      real M13 = n0 * c15[cell] * n0 + n0 * c14[cell] * n1 + n0 * c13[cell] * n2 + n1 * c56[cell] * n0 + n1 * c46[cell] * n1 + n1 * c36[cell] * n2 + n2 * c55[cell] * n0 + n2 * c35[cell] * n1 + n2 * c35[cell] * n2;
      real a[9];
      a[0] = M11;  
      a[1] = M12;  
      a[2] = M13;  
      a[3] = M12;  
      a[4] = M22;  
      a[5] = M23;  
      a[6] = M13;  
      a[7] = M23;  
      a[8] = M33;  
      Matrix33 A(a);
      saes.compute(A);
      auto eigenvalues = saes.eigenvalues();
      for(unsigned i = 0; i < 3; ++i) {
        maxEv = eigenvalues(i) > maxEv ? eigenvalues(i) : maxEv;
      }
    }
    pWaveVel[cell] = sqrt(maxEv / rho[cell]);
  }
  delete[] rho;
  delete[] c11;
  delete[] c12;
  delete[] c13;
  delete[] c14;
  delete[] c15;
  delete[] c16;
  delete[] c22;
  delete[] c23;
  delete[] c24;
  delete[] c25;
  delete[] c26;
  delete[] c33;
  delete[] c34;
  delete[] c35;
  delete[] c36;
  delete[] c44;
  delete[] c45;
  delete[] c46;
  delete[] c55;
  delete[] c56;
  delete[] c66;
#endif 
  double* timestep = new double[cells.size()];
  computeMaxTimesteps(mesh,pWaveVel, timestep);
  delete[] pWaveVel;  

  double localMinTimestep = *std::min_element(timestep, timestep + cells.size());
  double localMaxTimestep = *std::max_element(timestep, timestep + cells.size());
  double globalMinTimestep;
  double globalMaxTimestep;
#ifdef USE_MPI
  MPI_Allreduce(&localMinTimestep, &globalMinTimestep, 1, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm());
  MPI_Allreduce(&localMaxTimestep, &globalMaxTimestep, 1, MPI_DOUBLE, MPI_MAX, seissol::MPI::mpi.comm());
#else
  globalMinTimestep = localMinTimestep;
  globalMaxTimestep = localMaxTimestep;
#endif

  int* cluster = new int[cells.size()];
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    cluster[cell] = getCluster(timestep[cell], globalMinTimestep, m_rate);
  } 
  delete[] timestep;
  
  int totalNumberOfReductions = enforceMaximumDifference(mesh, cluster);

  delete[] m_vertexWeights;
  //m_ncon = 2;
  m_ncon = 1;
  m_vertexWeights = new int[cells.size() * m_ncon];
  int maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, m_rate);
  int drToCellRatio = 1;
  for (unsigned cell = 0; cell < cells.size(); ++cell) {    
    int dynamicRupture = 0;
    for (unsigned face = 0; face < 4; ++face) {
      dynamicRupture += ( getBoundaryCondition(boundaryCond, cell, face) == 3) ? 1 : 0;
    }
    
    m_vertexWeights[m_ncon * cell] = (1 + drToCellRatio*dynamicRupture) * ipow(m_rate, maxCluster - cluster[cell]);
    //m_vertexWeights[m_ncon * cell + 1] = (dynamicRupture > 0) ? 1 : 0;
  }

  delete[] cluster;

  logInfo(seissol::MPI::mpi.rank()) << "Computing LTS weights. Done. " << utils::nospace << '(' << totalNumberOfReductions << " reductions.)";
}

int seissol::initializers::time_stepping::LtsWeights::enforceMaximumDifference(PUML::TETPUML const& mesh, int* cluster) {
  int totalNumberOfReductions = 0;
  int globalNumberOfReductions;
  do {
    int localNumberOfReductions = enforceMaximumDifferenceLocal(mesh, cluster);

#ifdef USE_MPI
    MPI_Allreduce(&localNumberOfReductions, &globalNumberOfReductions, 1, MPI_INT, MPI_SUM, seissol::MPI::mpi.comm());    
#else
    globalNumberOfReductions = localNumberOfReductions;
#endif // USE_MPI
    totalNumberOfReductions += globalNumberOfReductions;
  } while (globalNumberOfReductions > 0);
  return totalNumberOfReductions;
}

int seissol::initializers::time_stepping::LtsWeights::enforceMaximumDifferenceLocal(PUML::TETPUML const& mesh, int* cluster, int maxDifference) {
  int numberOfReductions = 0;
  
  std::vector<PUML::TETPUML::cell_t> const& cells = mesh.cells();
  std::vector<PUML::TETPUML::face_t> const& faces = mesh.faces();
	int const* boundaryCond = mesh.cellData(1);

#ifdef USE_MPI
  std::unordered_map<int, std::vector<int>> rankToSharedFaces;
  std::unordered_map<int, int> localFaceIdToLocalCellId;
#endif // USE_MPI

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int timeCluster = cluster[cell];

		unsigned int faceids[4];
		PUML::Downward::faces(mesh, cells[cell], faceids);
    for (unsigned f = 0; f < 4; ++f) {
      int difference = maxDifference;
      int boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (boundary == 0 || boundary == 3 || boundary == 6) {
        // We treat MPI neighbours later
        auto const& face = faces[ faceids[f] ];
        if (!face.isShared()) {
          int cellIds[2];
          PUML::Upward::cells(mesh, face, cellIds);

          int neighbourCell = (cellIds[0] == static_cast<int>(cell)) ? cellIds[1] : cellIds[0];
          int otherTimeCluster = cluster[neighbourCell];
          
          if (boundary == 3) {
            difference = 0;
          }

          if (timeCluster > otherTimeCluster + difference) {
            timeCluster = otherTimeCluster + difference;
            ++numberOfReductions;
          }
        }
#ifdef USE_MPI
        else {
          rankToSharedFaces[ face.shared()[0] ].push_back(faceids[f]);
          localFaceIdToLocalCellId[ faceids[f] ] = cell;
        }
#endif // USE_MPI
      }
    }
    cluster[cell] = timeCluster;
  }

#ifdef USE_MPI
  FaceSorter faceSorter(faces);
  for (auto& sharedFaces: rankToSharedFaces) {
    std::sort(sharedFaces.second.begin(), sharedFaces.second.end(), faceSorter);
  }
  
  auto numExchanges = rankToSharedFaces.size();
  MPI_Request* requests = new MPI_Request[2*numExchanges];
  int** ghost = new int*[numExchanges];
  int** copy = new int*[numExchanges];
  auto exchange = rankToSharedFaces.begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.size();
    ghost[ex] = new int[ exchangeSize ];
    copy[ex]  = new int[ exchangeSize ];
    
    for (unsigned n = 0; n < exchangeSize; ++n) {
      copy[ex][n] = cluster[ localFaceIdToLocalCellId[ exchange->second[n] ] ];
    }
    MPI_Isend( copy[ex], exchangeSize, MPI_INT, exchange->first, 0, seissol::MPI::mpi.comm(), &requests[ex]);
    MPI_Irecv(ghost[ex], exchangeSize, MPI_INT, exchange->first, 0, seissol::MPI::mpi.comm(), &requests[numExchanges + ex]);
    ++exchange;
  }
  
  MPI_Waitall(2*numExchanges, requests, MPI_STATUSES_IGNORE);

  exchange = rankToSharedFaces.begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.size();
    for (unsigned n = 0; n < exchangeSize; ++n) {
      int difference = maxDifference;
      int otherTimeCluster = ghost[ex][n];
      
      int cellIds[2];
      PUML::Upward::cells(mesh, faces[ exchange->second[n] ], cellIds);
      int cell = (cellIds[0] >= 0) ? cellIds[0] : cellIds[1];

      unsigned int faceids[4];
      PUML::Downward::faces(mesh, cells[cell], faceids);
      unsigned f = 0;
      for (; f < 4 && static_cast<int>(faceids[f]) != exchange->second[n]; ++f);
      assert(f != 4);
      
      int boundary = getBoundaryCondition(boundaryCond, cell, f);
      if (boundary == 3) {
        difference = 0;
      }

      if (cluster[cell] > otherTimeCluster + difference) {
        cluster[cell] = otherTimeCluster + difference;
        ++numberOfReductions;
      }
    }
    ++exchange;
  }  
  
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    delete[] copy[ex];
    delete[] ghost[ex];
  }
  delete[] copy;
  delete[] ghost;
  delete[] requests;
#endif // USE_MPI

  return numberOfReductions;
}
