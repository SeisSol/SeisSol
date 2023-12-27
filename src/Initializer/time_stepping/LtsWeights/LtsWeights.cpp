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

#include <Initializer/time_stepping/GlobalTimestep.hpp>
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

double computeLocalCostOfClustering(const std::vector<int>& clusterIds,
                               const std::vector<int>& cellCosts,
                               unsigned int rate,
                               double wiggleFactor,
                               double minimalTimestep) {
  assert(clusterIds.size() == cellCosts.size());

  double cost = 0.0;
  for (auto i = 0U; i < clusterIds.size(); ++i) {
    const auto cluster = clusterIds[i];
    const auto cellCost = cellCosts[i];
    const double updateFactor = 1.0 / (std::pow(rate, cluster));
    cost += updateFactor * cellCost;
  }

  const auto minDtWithWiggle = minimalTimestep * wiggleFactor;
  return cost / minDtWithWiggle;
}

double computeGlobalCostOfClustering(const std::vector<int>& clusterIds,
                                     const std::vector<int>& cellCosts,
                                     unsigned int rate,
                                     double wiggleFactor,
                                     double minimalTimestep,
                                     MPI_Comm comm) {
  double cost =
      computeLocalCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, minimalTimestep);
#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &cost, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif // USE_MPI

  return cost;
}

std::vector<int> enforceMaxClusterId(const std::vector<int>& clusterIds, int maxClusterId) {
  auto newClusterIds = clusterIds;
  assert(maxClusterId >= 0);
  std::for_each(newClusterIds.begin(), newClusterIds.end(), [maxClusterId](auto& clusterId) {
    clusterId = std::min(clusterId, maxClusterId);
  });

  return newClusterIds;
}

// Merges clusters such that new cost is max oldCost * allowedPerformanceLossRatio
int computeMaxClusterIdAfterAutoMerge(const std::vector<int>& clusterIds,
                                      const std::vector<int>& cellCosts,
                                      unsigned int rate,
                                      double maximalAdmissibleCost,
                                      double wiggleFactor,
                                      double minimalTimestep) {
  int maxClusterId = *std::max_element(clusterIds.begin(), clusterIds.end());
#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, MPI::mpi.comm());
#endif

  // We only have one cluster for rate = 1 and thus cannot merge.
  if (rate == 1) {
    return maxClusterId;
  }

  // Iteratively merge clusters until we found the first number of clusters that has a cost that is too high
  for (auto curMaxClusterId = maxClusterId; curMaxClusterId >= 0; --curMaxClusterId) {
    const auto newClustering = enforceMaxClusterId(clusterIds, curMaxClusterId);
    const double cost = computeGlobalCostOfClustering(
        newClustering, cellCosts, rate, wiggleFactor, minimalTimestep, MPI::mpi.comm());
    if (cost > maximalAdmissibleCost) {
      // This is the first number of clusters that resulted in an inadmissible cost
      // Hence, it was admissible in the previous iteration
      return curMaxClusterId + 1;
    }
  }
  return 0;
}

LtsWeights::LtsWeights(const LtsWeightsConfig& config, seissol::SeisSol& seissolInstance)
    : m_velocityModel(config.velocityModel), m_rate(config.rate),
      m_vertexWeightElement(config.vertexWeightElement),
      m_vertexWeightDynamicRupture(config.vertexWeightDynamicRupture),
      m_vertexWeightFreeSurfaceWithGravity(config.vertexWeightFreeSurfaceWithGravity),
      seissolInstance(seissolInstance) { }

void LtsWeights::computeWeights(PUML::TETPUML const& mesh, double maximumAllowedTimeStep) {
  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Computing LTS weights.";

  // Note: Return value optimization is guaranteed while returning temp. objects in C++17
  m_mesh = &mesh;
  m_details = collectGlobalTimeStepDetails(maximumAllowedTimeStep);
  m_cellCosts = computeCostsPerTimestep();

  const auto* ltsParameters = seissolInstance.getMemoryManager().getLtsParameters();
  auto maxClusterIdToEnforce = ltsParameters->getMaxNumberOfClusters() - 1;
  if (ltsParameters->isWiggleFactorUsed() || ltsParameters->isAutoMergeUsed()) {
    auto autoMergeBaseline = ltsParameters->getAutoMergeCostBaseline();
    if (!(ltsParameters->isWiggleFactorUsed() && ltsParameters->isAutoMergeUsed())) {
      // Cost models only change things if both wiggle factor and auto merge are on.
      // In all other cases, choose the cheapest cost model.
      autoMergeBaseline = seissol::initializers::parameters::AutoMergeCostBaseline::MaxWiggleFactor;
    }

    ComputeWiggleFactorResult wiggleFactorResult{};
    if (autoMergeBaseline == seissol::initializers::parameters::AutoMergeCostBaseline::BestWiggleFactor) {
      // First compute wiggle factor without merging as baseline cost
      logInfo(rank) << "Using best wiggle factor as baseline cost for auto merging.";
      logInfo(rank) << "1. Compute best wiggle factor without merging clusters";
      const auto wiggleFactorResultBaseline = computeBestWiggleFactor(std::nullopt, false);
      // Compute wiggle factor a second time with merging and using the previous cost as baseline
      logInfo(rank) << "2. Compute best wiggle factor with merging clusters, using the previous cost estimate as baseline";
      const auto baselineCost = wiggleFactorResultBaseline.cost;
      wiggleFactorResult = computeBestWiggleFactor(baselineCost, ltsParameters->isAutoMergeUsed());
    } else {
      assert(autoMergeBaseline == AutoMergeCostBaseline::MaxWiggleFactor);
      wiggleFactorResult = computeBestWiggleFactor(std::nullopt, ltsParameters->isAutoMergeUsed());
    }

    wiggleFactor = wiggleFactorResult.wiggleFactor;
    if (ltsParameters->isAutoMergeUsed()) {
      maxClusterIdToEnforce =
          std::min(maxClusterIdToEnforce, wiggleFactorResult.maxClusterId);
    }
  } else {
    wiggleFactor = 1.0;
  }
  seissolInstance.wiggleFactorLts = wiggleFactor;

  m_clusterIds = computeClusterIds(wiggleFactor);

  m_ncon = evaluateNumberOfConstraints();
  auto finalNumberOfReductions = enforceMaximumDifference();

  logInfo(rank) << "Limiting number of clusters to" << maxClusterIdToEnforce + 1;
  m_clusterIds = enforceMaxClusterId(m_clusterIds, maxClusterIdToEnforce);

  int maxNumberOfClusters = *std::max_element(m_clusterIds.begin(), m_clusterIds.end()) + 1;
#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxNumberOfClusters, 1, MPI_INT, MPI_MAX, MPI::mpi.comm());
#endif
  seissolInstance.maxNumberOfClusters = maxNumberOfClusters;

  if (!m_vertexWeights.empty()) { m_vertexWeights.clear(); }
  m_vertexWeights.resize(m_clusterIds.size() * m_ncon);

  // calling virtual functions
  setVertexWeights();
  setAllowedImbalances();

  logInfo(rank) << "Computing LTS weights. Done. " << utils::nospace << '('
                                    << finalNumberOfReductions << " reductions.)";
}
LtsWeights::ComputeWiggleFactorResult
    LtsWeights::computeBestWiggleFactor(std::optional<double> baselineCost, bool isAutoMergeUsed) {
  const auto rank = seissol::MPI::mpi.rank();

  // Maps that keep track of number of clusters vs cost
  auto mapMaxClusterIdToLowestCost = std::map<int, double>{};
  auto maxMapClusterIdToBestWiggleFactor = std::map<int, double>{};

  const auto* ltsParameters = seissolInstance.getMemoryManager().getLtsParameters();
  const double minWiggleFactor = ltsParameters->getWiggleFactorMinimum();
  const double maxWiggleFactor = 1.0;

  const double stepSizeWiggleFactor = ltsParameters->getWiggleFactorStepsize();
  const int numberOfStepsWiggleFactor =
      std::ceil((maxWiggleFactor - minWiggleFactor) / stepSizeWiggleFactor) + 1;

  auto computeWiggleFactor = [minWiggleFactor, stepSizeWiggleFactor, maxWiggleFactor](auto ith) {
    return std::min(minWiggleFactor + ith * stepSizeWiggleFactor, maxWiggleFactor);
  };

  auto totalWiggleFactorReductions = 0u;

  if (baselineCost) {
    logInfo(rank) << "Baseline cost before cluster merging is" << *baselineCost;
  } else {
    // Compute baselineCost cost before wiggle factor and merging of clusters
    totalWiggleFactorReductions +=
        computeClusterIdsAndEnforceMaximumDifferenceCached(maxWiggleFactor);
    baselineCost = computeGlobalCostOfClustering(m_clusterIds,
                                                 m_cellCosts,
                                                 m_rate,
                                                 maxWiggleFactor,
                                                 m_details.globalMinTimeStep,
                                                 MPI::mpi.comm());
    logInfo(rank) << "Baseline cost, without wiggle factor and cluster merging is" << *baselineCost;
  }
  assert(baselineCost);

  const double maxAdmissibleCost =
      ltsParameters->getAllowedPerformanceLossRatioAutoMerge() * *baselineCost;

  if (isAutoMergeUsed) {
    logInfo(rank) << "Maximal admissible cost after cluster merging is" << maxAdmissibleCost;
  }

  for (int i = 0; i < numberOfStepsWiggleFactor; ++i) {
    const double curWiggleFactor = computeWiggleFactor(i);
    totalWiggleFactorReductions +=
        computeClusterIdsAndEnforceMaximumDifferenceCached(curWiggleFactor);

    // Note: Merging clusters does not invalidate invariance generated by enforceMaximumDifference()
    // This can be shown by enumerating all possible cases
    auto maxClusterIdToEnforce = ltsParameters->getMaxNumberOfClusters() - 1;
    if (isAutoMergeUsed) {
      const auto maxClusterIdAfterMerging =
          computeMaxClusterIdAfterAutoMerge(m_clusterIds,
                                            m_cellCosts,
                                            m_rate,
                                            maxAdmissibleCost,
                                            curWiggleFactor,
                                            m_details.globalMinTimeStep);
      maxClusterIdToEnforce = std::min(maxClusterIdAfterMerging, maxClusterIdToEnforce);
    }

    m_clusterIds = enforceMaxClusterId(m_clusterIds, maxClusterIdToEnforce);
    auto maxClusterId = *std::max_element(m_clusterIds.begin(), m_clusterIds.end());
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, MPI::mpi.comm());
#endif

    m_ncon = evaluateNumberOfConstraints();

    // Compute cost
    const double cost = computeGlobalCostOfClustering(m_clusterIds,
                                                      m_cellCosts,
                                                      m_rate,
                                                      curWiggleFactor,
                                                      m_details.globalMinTimeStep,
                                                      MPI::mpi.comm());

    if (auto it = mapMaxClusterIdToLowestCost.find(maxClusterId);
        it == mapMaxClusterIdToLowestCost.end() || cost <= it->second) {
      maxMapClusterIdToBestWiggleFactor[maxClusterId] = curWiggleFactor;
      mapMaxClusterIdToLowestCost[maxClusterId] = cost;
    }
  }

  // Find best wiggle factor after merging of clusters
  // We compare against cost of baselineCost.
  int minAdmissibleMaxClusterId = std::numeric_limits<int>::max();
  if (isAutoMergeUsed) {
    // When merging clusters, we want to find the minimum number of clusters with admissible performance.
    bool foundAdmissibleMerge = false;
    for (const auto& [noOfClusters, cost] : mapMaxClusterIdToLowestCost) {
      if (cost <= maxAdmissibleCost) {
        foundAdmissibleMerge = true;
        minAdmissibleMaxClusterId = std::min(minAdmissibleMaxClusterId, noOfClusters);
        logDebug(rank) << "Admissible. cluster:" << noOfClusters << ",cost" << cost
                       << "with wiggle factor" << maxMapClusterIdToBestWiggleFactor[noOfClusters];
      } else {
        logDebug(rank) << "Not admissible. cluster:" << noOfClusters << ",cost" << cost
                       << "with wiggle factor" << maxMapClusterIdToBestWiggleFactor[noOfClusters];
      }
    }
    if (!foundAdmissibleMerge) {
      logError() << "Found no admissible wiggle factor with cluster merge. Aborting.";
    }
  } else {
    // Otherwise choose one with the smallest cost
    minAdmissibleMaxClusterId =
        std::min_element(mapMaxClusterIdToLowestCost.begin(),
                         mapMaxClusterIdToLowestCost.end(),
                         [](const auto& a, const auto& b) { return a.second < b.second; })
            ->first;
  }

  logInfo(rank) << "Enforcing maximum difference when finding best wiggle factor took"
                << totalWiggleFactorReductions << "reductions.";

  const auto bestWiggleFactor = maxMapClusterIdToBestWiggleFactor[minAdmissibleMaxClusterId];
  const auto bestCostEstimate = mapMaxClusterIdToLowestCost[minAdmissibleMaxClusterId];
  logInfo(rank) << "The best wiggle factor is" << bestWiggleFactor << "with cost"
                << bestCostEstimate << "and" << minAdmissibleMaxClusterId + 1 << "time clusters";

  if (baselineCost > bestCostEstimate) {
    logInfo(rank) << "Cost decreased" << (*baselineCost - bestCostEstimate) / *baselineCost * 100
                  << "% with absolute cost decrease of" << *baselineCost - bestCostEstimate
                  << "compared to the baseline";
  } else {
    logInfo(rank) << "Cost increased" << (bestCostEstimate - *baselineCost) / *baselineCost * 100
                  << "% with absolute cost increase of" << bestCostEstimate - *baselineCost
                  << "compared to the baseline";
    logInfo(rank) << "Note: Cost increased due to cluster merging!";
  }

  return ComputeWiggleFactorResult{minAdmissibleMaxClusterId, bestWiggleFactor, bestCostEstimate};
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

seissol::initializers::GlobalTimestep LtsWeights::collectGlobalTimeStepDetails(double maximumAllowedTimeStep) {
  return seissol::initializers::computeTimesteps(1.0, maximumAllowedTimeStep, m_velocityModel, seissol::initializers::CellToVertexArray::fromPUML(*m_mesh), seissolInstance.getSeisSolParameters());
}

int LtsWeights::computeClusterIdsAndEnforceMaximumDifferenceCached(double curWiggleFactor) {
  int numberOfReductions = 0;
  auto lb = clusteringCache.lower_bound(curWiggleFactor);

  if (lb != clusteringCache.end() && !(clusteringCache.key_comp()(curWiggleFactor, lb->first))) {
    m_clusterIds = lb->second;
  } else {
    m_clusterIds = computeClusterIds(curWiggleFactor);
    const auto* ltsParameters = seissolInstance.getMemoryManager().getLtsParameters();
    if (ltsParameters->getWiggleFactorEnforceMaximumDifference()) {
      numberOfReductions = enforceMaximumDifference();
    }
    clusteringCache.insert(lb, std::make_pair(curWiggleFactor, m_clusterIds));
  }

  return numberOfReductions;
}

std::vector<int> LtsWeights::computeClusterIds(double curWiggleFactor) {
  const auto &cells = m_mesh->cells();
  std::vector<int> clusterIds(cells.size(), 0);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    clusterIds[cell] = getCluster(m_details.cellTimeStepWidths[cell],
                                  m_details.globalMinTimeStep, curWiggleFactor,
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
