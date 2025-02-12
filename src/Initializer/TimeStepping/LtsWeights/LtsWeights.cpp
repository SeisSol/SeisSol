// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "LtsWeights.h"

#include "Geometry/PUMLReader.h"
#include <Initializer/BasicTypedefs.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/Parameters/LtsParameters.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "PUML/Downward.h"
#include "PUML/PUML.h"
#include "PUML/Upward.h"

#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"
#include "generated_code/init.h"

namespace seissol::initializer::time_stepping {

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

  // Iteratively merge clusters until we found the first number of clusters that has a cost that is
  // too high
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
    : seissolInstance(seissolInstance), m_rate(config.rate),
      m_vertexWeightElement(config.vertexWeightElement),
      m_vertexWeightDynamicRupture(config.vertexWeightDynamicRupture),
      m_vertexWeightFreeSurfaceWithGravity(config.vertexWeightFreeSurfaceWithGravity),
      boundaryFormat(config.boundaryFormat) {}

void LtsWeights::computeWeights(PUML::TETPUML const& mesh) {
  logInfo() << "Computing LTS weights.";

  // Note: Return value optimization is guaranteed while returning temp. objects in C++17
  m_mesh = &mesh;
  m_details = collectGlobalTimeStepDetails();
  m_cellCosts = computeCostsPerTimestep();

  auto& ltsParameters = seissolInstance.getSeisSolParameters().timeStepping.lts;
  auto maxClusterIdToEnforce = ltsParameters.getMaxNumberOfClusters() - 1;

  prepareDifferenceEnforcement();

  if (ltsParameters.isWiggleFactorUsed() || ltsParameters.isAutoMergeUsed()) {
    auto autoMergeBaseline = ltsParameters.getAutoMergeCostBaseline();
    if (!(ltsParameters.isWiggleFactorUsed() && ltsParameters.isAutoMergeUsed())) {
      // Cost models only change things if both wiggle factor and auto merge are on.
      // In all other cases, choose the cheapest cost model.
      autoMergeBaseline = seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor;
    }

    ComputeWiggleFactorResult wiggleFactorResult{};
    if (autoMergeBaseline ==
        seissol::initializer::parameters::AutoMergeCostBaseline::BestWiggleFactor) {
      // First compute wiggle factor without merging as baseline cost
      logInfo() << "Using best wiggle factor as baseline cost for auto merging.";
      logInfo() << "1. Compute best wiggle factor without merging clusters";
      const auto wiggleFactorResultBaseline = computeBestWiggleFactor(std::nullopt, false);
      // Compute wiggle factor a second time with merging and using the previous cost as baseline
      logInfo() << "2. Compute best wiggle factor with merging clusters, using the previous cost "
                   "estimate as baseline";
      const auto baselineCost = wiggleFactorResultBaseline.cost;
      wiggleFactorResult = computeBestWiggleFactor(baselineCost, ltsParameters.isAutoMergeUsed());
    } else {
      assert(autoMergeBaseline ==
             seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor);
      wiggleFactorResult = computeBestWiggleFactor(std::nullopt, ltsParameters.isAutoMergeUsed());
    }

    wiggleFactor = wiggleFactorResult.wiggleFactor;
    if (ltsParameters.isAutoMergeUsed()) {
      maxClusterIdToEnforce = std::min(maxClusterIdToEnforce, wiggleFactorResult.maxClusterId);
    }
  } else {
    wiggleFactor = 1.0;
  }
  ltsParameters.setWiggleFactor(wiggleFactor);

  m_ncon = evaluateNumberOfConstraints();
  const auto finalNumberOfReductions =
      computeClusterIdsAndEnforceMaximumDifferenceCached(wiggleFactor);

  logInfo() << "Limiting number of clusters to" << maxClusterIdToEnforce + 1;
  m_clusterIds = enforceMaxClusterId(m_clusterIds, maxClusterIdToEnforce);

  int maxNumberOfClusters = *std::max_element(m_clusterIds.begin(), m_clusterIds.end()) + 1;
#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &maxNumberOfClusters, 1, MPI_INT, MPI_MAX, MPI::mpi.comm());
#endif
  ltsParameters.setMaxNumberOfClusters(maxNumberOfClusters);

  if (!m_vertexWeights.empty()) {
    m_vertexWeights.clear();
  }
  m_vertexWeights.resize(m_clusterIds.size() * m_ncon);

  // calling virtual functions
  setVertexWeights();
  setAllowedImbalances();

  logInfo() << "Computing LTS weights. Done. " << utils::nospace << '(' << finalNumberOfReductions
            << " reductions.)";
}
LtsWeights::ComputeWiggleFactorResult
    LtsWeights::computeBestWiggleFactor(std::optional<double> baselineCost, bool isAutoMergeUsed) {

  // Maps that keep track of number of clusters vs cost
  auto mapMaxClusterIdToLowestCost = std::map<int, double>{};
  auto maxMapClusterIdToBestWiggleFactor = std::map<int, double>{};

  const auto& ltsParameters = seissolInstance.getSeisSolParameters().timeStepping.lts;
  const double minWiggleFactor = ltsParameters.getWiggleFactorMinimum();
  const double maxWiggleFactor = 1.0;

  const double stepSizeWiggleFactor = ltsParameters.getWiggleFactorStepsize();
  const int numberOfStepsWiggleFactor =
      std::ceil((maxWiggleFactor - minWiggleFactor) / stepSizeWiggleFactor) + 1;

  auto computeWiggleFactor = [minWiggleFactor, stepSizeWiggleFactor, maxWiggleFactor](auto ith) {
    return std::min(minWiggleFactor + ith * stepSizeWiggleFactor, maxWiggleFactor);
  };

  auto totalWiggleFactorReductions = 0U;

  if (baselineCost) {
    logInfo() << "Baseline cost before cluster merging is" << *baselineCost;
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
    logInfo() << "Baseline cost, without wiggle factor and cluster merging is" << *baselineCost;
  }
  assert(baselineCost);

  const double maxAdmissibleCost =
      ltsParameters.getAllowedPerformanceLossRatioAutoMerge() * *baselineCost;

  if (isAutoMergeUsed) {
    logInfo() << "Maximal admissible cost after cluster merging is" << maxAdmissibleCost;
  }

  for (int i = 0; i < numberOfStepsWiggleFactor; ++i) {
    const double curWiggleFactor = computeWiggleFactor(i);
    totalWiggleFactorReductions +=
        computeClusterIdsAndEnforceMaximumDifferenceCached(curWiggleFactor);

    // Note: Merging clusters does not invalidate invariance generated by enforceMaximumDifference()
    // This can be shown by enumerating all possible cases
    auto maxClusterIdToEnforce = ltsParameters.getMaxNumberOfClusters() - 1;
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
    // When merging clusters, we want to find the minimum number of clusters with admissible
    // performance.
    bool foundAdmissibleMerge = false;
    for (const auto& [noOfClusters, cost] : mapMaxClusterIdToLowestCost) {
      if (cost <= maxAdmissibleCost) {
        foundAdmissibleMerge = true;
        minAdmissibleMaxClusterId = std::min(minAdmissibleMaxClusterId, noOfClusters);
        logDebug() << "Admissible. cluster:" << noOfClusters << ",cost" << cost
                   << "with wiggle factor" << maxMapClusterIdToBestWiggleFactor[noOfClusters];
      } else {
        logDebug() << "Not admissible. cluster:" << noOfClusters << ",cost" << cost
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

  logInfo() << "Enforcing maximum difference when finding best wiggle factor took"
            << totalWiggleFactorReductions << "reductions.";

  const auto bestWiggleFactor = maxMapClusterIdToBestWiggleFactor[minAdmissibleMaxClusterId];
  const auto bestCostEstimate = mapMaxClusterIdToLowestCost[minAdmissibleMaxClusterId];
  logInfo() << "The best wiggle factor is" << bestWiggleFactor << "with cost" << bestCostEstimate
            << "and" << minAdmissibleMaxClusterId + 1 << "time clusters";

  if (baselineCost > bestCostEstimate) {
    logInfo() << "Cost decreased" << (*baselineCost - bestCostEstimate) / *baselineCost * 100
              << "% with absolute cost decrease of" << *baselineCost - bestCostEstimate
              << "compared to the baseline";
  } else {
    logInfo() << "Cost increased" << (bestCostEstimate - *baselineCost) / *baselineCost * 100
              << "% with absolute cost increase of" << bestCostEstimate - *baselineCost
              << "compared to the baseline";
    logInfo() << "Note: Cost increased due to cluster merging!";
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

const std::vector<int>& LtsWeights::clusterIds() const { return m_clusterIds; }

const std::vector<double>& LtsWeights::timesteps() const { return m_details.cellTimeStepWidths; }

int LtsWeights::nWeightsPerVertex() const {
  assert(m_ncon != std::numeric_limits<int>::infinity() &&
         "num. constrains has not been initialized yet");
  return m_ncon;
}

int LtsWeights::getCluster(double timestep,
                           double globalMinTimestep,
                           double ltsWiggleFactor,
                           unsigned rate) {
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

FaceType LtsWeights::getBoundaryCondition(const void* boundaryCond, size_t cell, unsigned face) {
  int bcCurrentFace = seissol::geometry::decodeBoundary(boundaryCond, cell, face, boundaryFormat);
  if (bcCurrentFace > 64) {
    bcCurrentFace = 3;
  }
  return static_cast<FaceType>(bcCurrentFace);
}

int LtsWeights::ipow(int x, int y) {
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while (--y != 0) {
    result *= x;
  }
  return result;
}

seissol::initializer::GlobalTimestep LtsWeights::collectGlobalTimeStepDetails() {
  return seissol::initializer::computeTimesteps(
      seissol::initializer::CellToVertexArray::fromPUML(*m_mesh),
      seissolInstance.getSeisSolParameters());
}

int LtsWeights::computeClusterIdsAndEnforceMaximumDifferenceCached(double curWiggleFactor) {
  int numberOfReductions = 0;
  auto lb = clusteringCache.lower_bound(curWiggleFactor);

  if (lb != clusteringCache.end() && !(clusteringCache.key_comp()(curWiggleFactor, lb->first))) {
    m_clusterIds = lb->second;
  } else {
    // re-use best computed maxdiff enforcement available
    // reason that works: cf. Lukas' proof for cluster merging not violating maximum difference
    // we may generalize due to the fact that min(a, min(b,c)) = min(min(a,b), c) = min(min(a,c),
    // b), essentially establishing a partial ordering of clusterings, where A >= B iff
    // cluster(A[i]) >= cluster(B[i]) for all cells i. Thus: walking through the wiggle factors from
    // lower to higher will save a lot of reductions

    int cellchanges = 0;
    if (lb != clusteringCache.end()) {
      // use the cache
      const auto newClusterIds = computeClusterIds(curWiggleFactor);
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : cellchanges)
#endif
      for (unsigned cell = 0; cell < m_mesh->cells().size(); ++cell) {
        if (lb->second[cell] > newClusterIds[cell]) {
          ++cellchanges;
        }
        m_clusterIds[cell] = std::min(lb->second[cell], newClusterIds[cell]);
      }
    } else {
      m_clusterIds = computeClusterIds(curWiggleFactor);
      cellchanges = m_mesh->cells().size();
    }
    const auto& ltsParameters = seissolInstance.getSeisSolParameters().timeStepping.lts;
    if (ltsParameters.getWiggleFactorEnforceMaximumDifference()) {
#ifdef USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &cellchanges, 1, MPI_INT, MPI_SUM, seissol::MPI::mpi.comm());
#endif
      if (cellchanges > 0) {
        numberOfReductions = enforceMaximumDifference();
      }
    }
    clusteringCache[curWiggleFactor] = m_clusterIds;
  }

  return numberOfReductions;
}

std::vector<int> LtsWeights::computeClusterIds(double curWiggleFactor) {
  const auto& cells = m_mesh->cells();
  std::vector<int> clusterIds(cells.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    clusterIds[cell] = getCluster(
        m_details.cellTimeStepWidths[cell], m_details.globalMinTimeStep, curWiggleFactor, m_rate);
  }
  return clusterIds;
}

std::vector<int> LtsWeights::computeCostsPerTimestep() {
  const auto& cells = m_mesh->cells();

  std::vector<int> cellCosts(cells.size());
  const void* boundaryCond = m_mesh->cellData(1);
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int dynamicRupture = 0;
    int freeSurfaceWithGravity = 0;

    unsigned int faceids[4];
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);

    for (unsigned face = 0; face < 4; ++face) {
      const auto faceType = getBoundaryCondition(boundaryCond, cell, face);
      dynamicRupture += (faceType == FaceType::DynamicRupture) ? 1 : 0;
      freeSurfaceWithGravity += (faceType == FaceType::FreeSurfaceGravity) ? 1 : 0;
    }

    const int costDynamicRupture = m_vertexWeightDynamicRupture * dynamicRupture;
    const int costDisplacement = m_vertexWeightFreeSurfaceWithGravity * freeSurfaceWithGravity;
    cellCosts[cell] = m_vertexWeightElement + costDynamicRupture + costDisplacement;
  }
  return cellCosts;
}

int LtsWeights::enforceMaximumDifference() {
  int totalNumberOfReductions = 0;
  int globalNumberOfReductions = 0;
  do {
    int localNumberOfReductions = enforceMaximumDifferenceLocal();

#ifdef USE_MPI
    MPI_Allreduce(&localNumberOfReductions,
                  &globalNumberOfReductions,
                  1,
                  MPI_INT,
                  MPI_SUM,
                  seissol::MPI::mpi.comm());
#else
    globalNumberOfReductions = localNumberOfReductions;
#endif // USE_MPI
    totalNumberOfReductions += globalNumberOfReductions;
  } while (globalNumberOfReductions > 0);
  return totalNumberOfReductions;
}

void LtsWeights::prepareDifferenceEnforcement() {
#ifdef USE_MPI
  const auto& cells = m_mesh->cells();
  const auto& faces = m_mesh->faces();
  const void* boundaryCond = m_mesh->cellData(1);

  std::unordered_map<int, std::vector<int>> rankToSharedFacesPre;
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    unsigned int faceids[4]{};
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);
    for (unsigned f = 0; f < 4; ++f) {
      const auto boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbours later
        const auto& face = faces.at(faceids[f]);
        if (face.isShared()) {
          rankToSharedFacesPre[face.shared()[0]].push_back(faceids[f]);
          localFaceIdToLocalCellId[faceids[f]] = cell;
        }
      }
    }
  }

  for (auto& sharedFaces : rankToSharedFacesPre) {
    std::sort(sharedFaces.second.begin(),
              sharedFaces.second.end(),
              [&](unsigned int a, unsigned int b) { return faces[a].gid() < faces[b].gid(); });
  }

  rankToSharedFaces =
      decltype(rankToSharedFaces)(rankToSharedFacesPre.begin(), rankToSharedFacesPre.end());
#endif // USE_MPI
}

int LtsWeights::enforceMaximumDifferenceLocal(int maxDifference) {
  int numberOfReductions = 0;

  const auto& cells = m_mesh->cells();
  const auto& faces = m_mesh->faces();
  const void* boundaryCond = m_mesh->cellData(1);

  const auto cellCount = cells.size();

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numberOfReductions)
#endif
  for (unsigned cell = 0; cell < cells.size(); ++cell) {
    int timeCluster = m_clusterIds[cell];

    unsigned int faceids[4]{};
    PUML::Downward::faces(*m_mesh, cells[cell], faceids);
    for (unsigned f = 0; f < 4; ++f) {
      int difference = maxDifference;
      const auto boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbours later
        const auto& face = faces.at(faceids[f]);
        if (!face.isShared()) {
          int cellIds[2];
          PUML::Upward::cells(*m_mesh, face, cellIds);

          const int neighborCell = (cellIds[0] == static_cast<int>(cell)) ? cellIds[1] : cellIds[0];
          const int otherTimeCluster = m_clusterIds[neighborCell];

          if (boundary == FaceType::DynamicRupture) {
            difference = 0;
          }

          if (timeCluster > otherTimeCluster + difference) {
            timeCluster = otherTimeCluster + difference;
            ++numberOfReductions;
          }
        }
      }
    }
    m_clusterIds[cell] = timeCluster;
  }

#ifdef USE_MPI
  const auto numExchanges = rankToSharedFaces.size();
  std::vector<MPI_Request> requests(2 * numExchanges);
  std::vector<std::vector<int>> ghost(numExchanges);
  std::vector<std::vector<int>> copy(numExchanges);

  for (std::size_t ex = 0; ex < numExchanges; ++ex) {
    const auto& exchange = rankToSharedFaces[ex];
    const auto exchangeSize = exchange.second.size();
    ghost[ex].resize(exchangeSize);
    copy[ex].resize(exchangeSize);

    for (std::size_t n = 0; n < exchangeSize; ++n) {
      copy[ex][n] = m_clusterIds[localFaceIdToLocalCellId[exchange.second[n]]];
    }
    MPI_Isend(copy[ex].data(),
              exchangeSize,
              MPI_INT,
              exchange.first,
              0,
              seissol::MPI::mpi.comm(),
              &requests[ex]);
    MPI_Irecv(ghost[ex].data(),
              exchangeSize,
              MPI_INT,
              exchange.first,
              0,
              seissol::MPI::mpi.comm(),
              &requests[numExchanges + ex]);
  }

  MPI_Waitall(2 * numExchanges, requests.data(), MPI_STATUSES_IGNORE);

  auto* idData = m_clusterIds.data();
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numberOfReductions) reduction(min : idData[0 : cellCount])
#endif
  for (std::size_t ex = 0; ex < numExchanges; ++ex) {
    const auto& exchange = rankToSharedFaces[ex];
    const auto exchangeSize = exchange.second.size();
    for (std::size_t n = 0; n < exchangeSize; ++n) {
      int difference = maxDifference;
      const int otherTimeCluster = ghost[ex][n];

      int cellIds[2];
      PUML::Upward::cells(*m_mesh, faces[exchange.second[n]], cellIds);
      const int cell = (cellIds[0] >= 0) ? cellIds[0] : cellIds[1];

      unsigned int faceids[4];
      PUML::Downward::faces(*m_mesh, cells[cell], faceids);
      unsigned f = 0;
      for (; f < 4 && static_cast<int>(faceids[f]) != exchange.second[n]; ++f) {
      }
      assert(f != 4);

      const auto boundary = getBoundaryCondition(boundaryCond, cell, f);
      if (boundary == FaceType::DynamicRupture) {
        difference = 0;
      }

      if (m_clusterIds[cell] > otherTimeCluster + difference) {
        m_clusterIds[cell] = otherTimeCluster + difference;
        ++numberOfReductions;
      }
    }
  }

#endif // USE_MPI

  return numberOfReductions;
}
} // namespace seissol::initializer::time_stepping
