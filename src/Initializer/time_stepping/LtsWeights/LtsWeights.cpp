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

#include <Initializer/ParameterDB.h>
#include <Parallel/MPI.h>

#include <generated_code/init.h>

#include "SeisSol.h"

namespace seissol::initializers::time_stepping {

class FaceSorter {
private:
  std::vector<PUML::TETPUML::face_t> const &m_faces;

public:
  FaceSorter(std::vector<PUML::TETPUML::face_t> const &faces) : m_faces(faces) {}

  bool operator()(unsigned int a, unsigned int b) const {
    return m_faces[a].gid() < m_faces[b].gid();
  }
};

void LtsWeights::computeWeights(PUML::TETPUML const &mesh, double maximumAllowedTimeStep) {
  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Computing LTS weights.";

  // Note: Return value optimization is guaranteed while returning temp. objects in C++17
  m_mesh = &mesh;
  m_details = collectGlobalTimeStepDetails(maximumAllowedTimeStep);
  m_cellCosts = computeCostsPerTimestep();

  const double stepSizeWiggleFactor = 0.01;
  const double minWiggleFactor = 1.0 / m_rate + stepSizeWiggleFactor;
  const double maxWiggleFactor = 1.0;
  const int numberOfStepsWiggleFactor = std::ceil((maxWiggleFactor - minWiggleFactor)/ stepSizeWiggleFactor) + 1;

  auto computeWiggleFactor = [minWiggleFactor, stepSizeWiggleFactor, maxWiggleFactor](auto ith) {
    return std::min(minWiggleFactor + ith * stepSizeWiggleFactor, maxWiggleFactor);
  };

  auto costEstimates = std::vector<double>(numberOfStepsWiggleFactor, std::numeric_limits<double>::max());

  auto totalWiggleFactorReductions = 0u;
  for (int i = 0; i <= numberOfStepsWiggleFactor; ++i) {
    const double curWiggleFactor = computeWiggleFactor(i);
    m_clusterIds = computeClusterIds(curWiggleFactor);
    m_ncon = evaluateNumberOfConstraints();
    totalWiggleFactorReductions += enforceMaximumDifference();

    // Compute cost
    costEstimates[i] = 0.0;
    for (auto j=0U; j < m_clusterIds.size(); ++j) {
      auto cluster = m_clusterIds[j];
      const double updateFactor = 1.0/(std::pow(2, cluster) * curWiggleFactor * m_details.globalMinTimeStep);

      costEstimates[i] += updateFactor * m_cellCosts[j];
    }
  }
#ifdef USE_MPI
  MPI_Allreduce(
      MPI_IN_PLACE,
      costEstimates.data(),
      static_cast<int>(costEstimates.size()),
      MPI_DOUBLE,
      MPI_SUM,
      seissol::MPI::mpi.comm()
      );
  MPI_Barrier(seissol::MPI::mpi.comm());
#endif

  auto maxWiggleFactorCostEstimate = costEstimates[costEstimates.size() - 1];

  double bestCostEstimate = std::numeric_limits<double>::max();
  double bestWiggleFactor = maxWiggleFactor;
  for (auto i = 0u; i < costEstimates.size(); ++i) {
    auto curCost = costEstimates[i];
    logDebug(rank) << "A wiggle factor of "
        << computeWiggleFactor(i)
        << "leads to a cost of"
        << curCost
        << "which is a speedup  of"
        << (maxWiggleFactorCostEstimate / curCost) * 100 - 100
        << "%.";
    // Note: Higher wiggle factor with same cost is better!
    if (curCost <= bestCostEstimate) {
      bestCostEstimate = curCost;
      bestWiggleFactor = computeWiggleFactor(i);
    }
  }

  logInfo(rank) << "Best wiggle factor" << bestWiggleFactor << "with cost" << bestCostEstimate;
  logInfo(rank) << "Finding best factor took" << totalWiggleFactorReductions << "reductions.";

  logInfo(rank) << "Speedup of" << (maxWiggleFactorCostEstimate / bestCostEstimate) * 100 - 100
      << "% with absolute cost difference" << maxWiggleFactorCostEstimate - bestCostEstimate
      << "compared to the default wiggle factor of"
      << maxWiggleFactor;

  seissol::SeisSol::main.wiggleFactorLts = bestWiggleFactor;
  wiggleFactor = bestWiggleFactor;

  m_clusterIds = computeClusterIds(wiggleFactor);

  m_ncon = evaluateNumberOfConstraints();
  auto finalNumberOfReductions = enforceMaximumDifference();

  if (!m_vertexWeights.empty()) { m_vertexWeights.clear(); }
  m_vertexWeights.resize(m_clusterIds.size() * m_ncon);

  // calling virtual functions
  setVertexWeights();
  setAllowedImbalances();

  logInfo(rank) << "Computing LTS weights. Done. " << utils::nospace << '('
                                    << finalNumberOfReductions << " reductions.)";
}

const int* LtsWeights::vertexWeights() const {
  assert(!m_vertexWeights.empty() && "vertex weights are not initialized");
  return m_vertexWeights.data();
}

const double* LtsWeights::imbalances() const {
  assert(!m_imbalances.empty() && "weight imbalances are not initialized");
  return m_imbalances.data();
}

int LtsWeights::nWeightsPerVertex() const {
  assert(m_ncon != std::numeric_limits<int>::infinity() && "num. constrains has not been initialized yet");
  return m_ncon;
}

void LtsWeights::computeMaxTimesteps(std::vector<double> const &pWaveVel,
                                     std::vector<double> &timeSteps, double maximumAllowedTimeStep) {
  std::vector<PUML::TETPUML::cell_t> const &cells = m_mesh->cells();
  std::vector<PUML::TETPUML::vertex_t> const &vertices = m_mesh->vertices();

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    // Compute insphere radius
    Eigen::Vector3d barycentre(0., 0., 0.);
    Eigen::Vector3d x[4];
    unsigned vertLids[4];
    PUML::Downward::vertices(*m_mesh, cells[cell], vertLids);
    for (unsigned vtx = 0; vtx < 4; ++vtx) {
      for (unsigned d = 0; d < 3; ++d) {
        x[vtx](d) = vertices[vertLids[vtx]].coordinate()[d];
      }
    }
    Eigen::Matrix4d A;
    A << x[0](0), x[0](1), x[0](2), 1.0,
        x[1](0), x[1](1), x[1](2), 1.0,
        x[2](0), x[2](1), x[2](2), 1.0,
        x[3](0), x[3](1), x[3](2), 1.0;

    double alpha = A.determinant();
    double Nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
    double Nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
    double Nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
    double Nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
    double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);

    // Compute maximum timestep (CFL=1)
    timeSteps[cell] = std::fmin(maximumAllowedTimeStep, 2.0 * insphere / (pWaveVel[cell] * (2 * CONVERGENCE_ORDER - 1)));
  }
}

int LtsWeights::getCluster(double timestep, double globalMinTimestep, double ltsWiggleFactor, unsigned rate) {
  if (rate == 1) {
    return 0;
  }

  double upper = ltsWiggleFactor * rate * globalMinTimestep;

  int cluster = 0;
  while (upper <= timestep) {
    upper *= rate;
    ++cluster;
  }
  return cluster;
}

int LtsWeights::getBoundaryCondition(int const *boundaryCond, unsigned cell, unsigned face) {
  int bcCurrentFace = ((boundaryCond[cell] >> (face * 8)) & 0xFF);
  if (bcCurrentFace > 64) {
    bcCurrentFace = 3;
  }
  return bcCurrentFace;
}

int LtsWeights::ipow(int x, int y) {
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while (--y) {
    result *= x;
  }
  return result;
}

LtsWeights::GlobalTimeStepDetails LtsWeights::collectGlobalTimeStepDetails(double maximumAllowedTimeStep) {

  const auto &cells = m_mesh->cells();
  std::vector<double> pWaveVel;
  pWaveVel.resize(cells.size());

  // TODO(Sebastian) Use averaged generator here as well.
  seissol::initializers::ElementBarycentreGeneratorPUML queryGen(*m_mesh);
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
  parameterDB.evaluateModel(m_velocityModel, &queryGen);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    pWaveVel[cell] = materials[cell].getMaxWaveSpeed();
  }

  GlobalTimeStepDetails details{};
  details.timeSteps.resize(cells.size());
  computeMaxTimesteps(pWaveVel, details.timeSteps, maximumAllowedTimeStep);

  double localMinTimestep = *std::min_element(details.timeSteps.begin(), details.timeSteps.end());
  double localMaxTimestep = *std::max_element(details.timeSteps.begin(), details.timeSteps.end());

#ifdef USE_MPI
  MPI_Allreduce(&localMinTimestep, &details.globalMinTimeStep, 1, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm());
  MPI_Allreduce(&localMaxTimestep, &details.globalMaxTimeStep, 1, MPI_DOUBLE, MPI_MAX, seissol::MPI::mpi.comm());
#else
  details.globalMinTimeStep = localMinTimestep;
  details.globalMaxTimeStep = localMaxTimestep;
#endif
  return details;
}

std::vector<int> LtsWeights::computeClusterIds(double wiggleFactor) {
  const auto &cells = m_mesh->cells();
  std::vector<int> clusterIds(cells.size(), 0);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    clusterIds[cell] = getCluster(m_details.timeSteps[cell],
                                  m_details.globalMinTimeStep,
                                  wiggleFactor,
                                  m_rate);
  }
  return clusterIds;
}

std::vector<int> LtsWeights::computeCostsPerTimestep() {
  const auto &cells = m_mesh->cells();

  std::vector<int> cellCosts(cells.size());
  int const *boundaryCond = m_mesh->cellData(1);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int dynamicRupture = 0;
    int freeSurfaceWithGravity = 0;

    unsigned int faceids[4];
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);

    for (unsigned face = 0; face < 4; ++face) {
      const auto faceType = static_cast<FaceType>(getBoundaryCondition(boundaryCond, cell, face));
      dynamicRupture += (faceType == FaceType::dynamicRupture) ? 1 : 0;
      freeSurfaceWithGravity += (faceType == FaceType::freeSurfaceGravity) ? 1 : 0;
    }

    const int costDynamicRupture = m_vertexWeightDynamicRupture * dynamicRupture;
    const int costDisplacement = m_vertexWeightFreeSurfaceWithGravity * freeSurfaceWithGravity;
    cellCosts[cell] = m_vertexWeightElement + costDynamicRupture + costDisplacement;
  }
  return cellCosts;
}

int LtsWeights::enforceMaximumDifference() {
  int totalNumberOfReductions = 0;
  int globalNumberOfReductions;
  do {
    int localNumberOfReductions = enforceMaximumDifferenceLocal();

#ifdef USE_MPI
    MPI_Allreduce(&localNumberOfReductions, &globalNumberOfReductions, 1, MPI_INT, MPI_SUM, seissol::MPI::mpi.comm());
#else
    globalNumberOfReductions = localNumberOfReductions;
#endif // USE_MPI
    totalNumberOfReductions += globalNumberOfReductions;
  } while (globalNumberOfReductions > 0);
  return totalNumberOfReductions;
}

int LtsWeights::enforceMaximumDifferenceLocal(int maxDifference) {
  int numberOfReductions = 0;

  std::vector<PUML::TETPUML::cell_t> const &cells = m_mesh->cells();
  std::vector<PUML::TETPUML::face_t> const &faces = m_mesh->faces();
  int const *boundaryCond = m_mesh->cellData(1);

#ifdef USE_MPI
  std::unordered_map<int, std::vector<int>> rankToSharedFaces;
  std::unordered_map<int, int> localFaceIdToLocalCellId;
#endif // USE_MPI

  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int timeCluster = m_clusterIds[cell];

    unsigned int faceids[4];
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);
    for (unsigned f = 0; f < 4; ++f) {
      int difference = maxDifference;
      int boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (boundary == 0 || boundary == 3 || boundary == 6) {
        // We treat MPI neighbours later
        auto const &face = faces[faceids[f]];
        if (!face.isShared()) {
          int cellIds[2];
          PUML::Upward::cells(*m_mesh, face, cellIds);

          int neighbourCell = (cellIds[0] == static_cast<int>(cell)) ? cellIds[1] : cellIds[0];
          int otherTimeCluster = m_clusterIds[neighbourCell];

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
          rankToSharedFaces[face.shared()[0]].push_back(faceids[f]);
          localFaceIdToLocalCellId[faceids[f]] = cell;
        }
#endif // USE_MPI
      }
    }
    m_clusterIds[cell] = timeCluster;
  }

#ifdef USE_MPI
  FaceSorter faceSorter(faces);
  for (auto &sharedFaces: rankToSharedFaces) {
    std::sort(sharedFaces.second.begin(), sharedFaces.second.end(), faceSorter);
  }

  auto numExchanges = rankToSharedFaces.size();
  std::vector<MPI_Request> requests(2 * numExchanges);
  std::vector<std::vector<int>> ghost(numExchanges);
  std::vector<std::vector<int>> copy(numExchanges);

  auto exchange = rankToSharedFaces.begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.size();
    ghost[ex].resize(exchangeSize);
    copy[ex].resize(exchangeSize);

    for (unsigned n = 0; n < exchangeSize; ++n) {
      copy[ex][n] = m_clusterIds[localFaceIdToLocalCellId[exchange->second[n]]];
    }
    MPI_Isend(copy[ex].data(), exchangeSize, MPI_INT, exchange->first, 0, seissol::MPI::mpi.comm(), &requests[ex]);
    MPI_Irecv(ghost[ex].data(), exchangeSize, MPI_INT, exchange->first, 0, seissol::MPI::mpi.comm(),
              &requests[numExchanges + ex]);
    ++exchange;
  }

  MPI_Waitall(2 * numExchanges, requests.data(), MPI_STATUSES_IGNORE);

  exchange = rankToSharedFaces.begin();
  for (unsigned ex = 0; ex < numExchanges; ++ex) {
    auto exchangeSize = exchange->second.size();
    for (unsigned n = 0; n < exchangeSize; ++n) {
      int difference = maxDifference;
      int otherTimeCluster = ghost[ex][n];

      int cellIds[2];
      PUML::Upward::cells(*m_mesh, faces[exchange->second[n]], cellIds);
      int cell = (cellIds[0] >= 0) ? cellIds[0] : cellIds[1];

      unsigned int faceids[4];
      PUML::Downward::faces(*m_mesh, cells[cell], faceids);
      unsigned f = 0;
      for (; f < 4 && static_cast<int>(faceids[f]) != exchange->second[n]; ++f);
      assert(f != 4);

      int boundary = getBoundaryCondition(boundaryCond, cell, f);
      if (boundary == 3) {
        difference = 0;
      }

      if (m_clusterIds[cell] > otherTimeCluster + difference) {
        m_clusterIds[cell] = otherTimeCluster + difference;
        ++numberOfReductions;
      }
    }
    ++exchange;
  }

#endif // USE_MPI

  return numberOfReductions;
}
} // namespace seissol::initializers::time_stepping
