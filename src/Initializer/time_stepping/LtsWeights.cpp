/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017 - 2020, SeisSol Group
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
#include <Initializer/typedefs.hpp>

#include <PUML/PUML.h>
#include <PUML/Downward.h>
#include <PUML/Upward.h>
#include "LtsWeights.h"

#include <Eigen/Dense>

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
                                                                            std::vector<double> const& pWaveVel,
                                                                            std::vector<double>& timestep ) {
  std::vector<PUML::TETPUML::cell_t> const& cells = mesh.cells();
  std::vector<PUML::TETPUML::vertex_t> const& vertices = mesh.vertices();

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    // Compute insphere radius
    Eigen::Vector3d barycentre(0.,0.,0.);
    Eigen::Vector3d x[4];
    unsigned vertLids[4];
    PUML::Downward::vertices(mesh, cells[cell], vertLids);
    for (unsigned vtx = 0; vtx < 4; ++vtx) {
      for (unsigned d = 0; d < 3; ++d) {
        x[vtx](d) = vertices[ vertLids[vtx] ].coordinate()[d];
      }
    }
    Eigen::Matrix4d A;
    A << x[0](0), x[0](1), x[0](2), 1.0,
         x[1](0), x[1](1), x[1](2), 1.0,
         x[2](0), x[2](1), x[2](2), 1.0,
         x[3](0), x[3](1), x[3](2), 1.0;

    double alpha = A.determinant();
    double Nabc = ( (x[1]-x[0]).cross(x[2]-x[0]) ).norm();
    double Nabd = ( (x[1]-x[0]).cross(x[3]-x[0]) ).norm();
    double Nacd = ( (x[2]-x[0]).cross(x[3]-x[0]) ).norm();
    double Nbcd = ( (x[2]-x[1]).cross(x[3]-x[1]) ).norm();
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
  if (bcCurrentFace > 64) {
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

  const auto& cells = mesh.cells();
  const auto& faces = mesh.faces();
  int const* boundaryCond = mesh.cellData(1);

  std::vector<double> pWaveVel;
  pWaveVel.resize(cells.size());
  
  seissol::initializers::ElementBarycentreGeneratorPUML queryGen(mesh);  
  //up to now we only distinguish between anisotropic elastic any other isotropic material
#ifdef USE_ANISOTROPIC
  std::vector<seissol::model::AnisotropicMaterial> materials(cells.size());
  seissol::initializers::MaterialParameterDB<seissol::model::AnisotropicMaterial> parameterDB;
#elif defined(USE_POROELASTIC)
  std::vector<seissol::model::PoroElasticMaterial> materials(cells.size());
  seissol::initializers::MaterialParameterDB<seissol::model::PoroElasticMaterial> parameterDB;
#else
  std::vector<seissol::model::ElasticMaterial> materials(cells.size());
  seissol::initializers::MaterialParameterDB<seissol::model::ElasticMaterial> parameterDB;
#endif 
  parameterDB.setMaterialVector(&materials);
  parameterDB.evaluateModel(m_velocityModel, queryGen);
  for(unsigned cell = 0; cell < cells.size(); ++cell) {
    pWaveVel[cell] = materials[cell].getMaxWaveSpeed();
  }
  std::vector<double> timestep;
  timestep.resize(cells.size());
  computeMaxTimesteps(mesh, pWaveVel, timestep);

  double localMinTimestep = *std::min_element(timestep.begin(), timestep.end());
  double localMaxTimestep = *std::max_element(timestep.begin(), timestep.end());
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
  
  int totalNumberOfReductions = enforceMaximumDifference(mesh, cluster);

  delete[] m_vertexWeights;
  //m_ncon = 2;
  m_ncon = 1;
  m_vertexWeights = new int[cells.size() * m_ncon];
  int maxCluster = getCluster(globalMaxTimestep, globalMinTimestep, m_rate);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int dynamicRupture = 0;
    int freeSurface = 0;

    unsigned int faceids[4];
    PUML::Downward::faces(mesh, cells[cell], faceids);

    for (unsigned face = 0; face < 4; ++face) {
      const auto faceType = static_cast<FaceType>(getBoundaryCondition(boundaryCond, cell, face));
      dynamicRupture += ( getBoundaryCondition(boundaryCond, cell, face) == 3) ? 1 : 0;

      int cellIds[2];
      PUML::Upward::cells(mesh, faces[faceids[face]], cellIds);

      if (faceType == FaceType::freeSurfaceGravity) {
        freeSurface++;
      }
    }

    const int costDynamicRupture = vertexWeightDynamicRupture * dynamicRupture;
    const int costDisplacement = vertexWeightFreeSurfaceWithGravity * freeSurface;
    const int costPerTimestep = vertexWeightElement + costDynamicRupture + costDisplacement;
    m_vertexWeights[m_ncon * cell] = costPerTimestep * ipow(m_rate, maxCluster - cluster[cell]);
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
