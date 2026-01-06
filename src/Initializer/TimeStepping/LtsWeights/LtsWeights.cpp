// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "LtsWeights.h"

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"

#include <PUML/Downward.h>
#include <PUML/Upward.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <mpi.h>
#include <optional>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::time_stepping {

double computeLocalCostOfClustering(const std::vector<int>& clusterIds,
                                    const std::vector<int>& cellCosts,
                                    const std::vector<uint64_t>& rate,
                                    double wiggleFactor,
                                    double minimalTimestep) {
  assert(clusterIds.size() == cellCosts.size());

  double cost = 0.0;
  for (auto i = 0U; i < clusterIds.size(); ++i) {
    const auto cluster = clusterIds[i];
    const auto cellCost = cellCosts[i];
    const auto invUpdateFactor = ratepow(rate, 0, cluster);
    cost += cellCost / static_cast<double>(invUpdateFactor);
  }

  const auto minDtWithWiggle = minimalTimestep * wiggleFactor;
  return cost / minDtWithWiggle;
}

double computeGlobalCostOfClustering(const std::vector<int>& clusterIds,
                                     const std::vector<int>& cellCosts,
                                     const std::vector<uint64_t>& rate,
                                     double wiggleFactor,
                                     double minimalTimestep,
                                     MPI_Comm comm) {
  double cost =
      computeLocalCostOfClustering(clusterIds, cellCosts, rate, wiggleFactor, minimalTimestep);
  MPI_Allreduce(MPI_IN_PLACE, &cost, 1, MPI_DOUBLE, MPI_SUM, comm);

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
                                      const std::vector<uint64_t>& rate,
                                      double maximalAdmissibleCost,
                                      double wiggleFactor,
                                      double minimalTimestep) {
  int maxClusterId = *std::max_element(clusterIds.begin(), clusterIds.end());
  MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, Mpi::mpi.comm());

  // Iteratively merge clusters until we found the first number of clusters that has a cost that is
  // too high
  for (auto curMaxClusterId = maxClusterId; curMaxClusterId >= 0; --curMaxClusterId) {
    const auto newClustering = enforceMaxClusterId(clusterIds, curMaxClusterId);
    const double cost = computeGlobalCostOfClustering(
        newClustering, cellCosts, rate, wiggleFactor, minimalTimestep, Mpi::mpi.comm());
    if (cost > maximalAdmissibleCost) {
      // This is the first number of clusters that resulted in an inadmissible cost
      // Hence, it was admissible in the previous iteration
      return std::min(maxClusterId, curMaxClusterId + 1);
    }
  }
  return 0;
}

LtsWeights::LtsWeights(const LtsWeightsConfig& config, seissol::SeisSol& seissolInstance)
    : seissolInstance_(seissolInstance), rate_(config.rate),
      vertexWeightElement_(config.vertexWeightElement),
      vertexWeightDynamicRupture_(config.vertexWeightDynamicRupture),
      vertexWeightFreeSurfaceWithGravity_(config.vertexWeightFreeSurfaceWithGravity),
      boundaryFormat_(config.boundaryFormat) {}

void LtsWeights::computeWeights(const seissol::geometry::PumlMesh& meshTopology,
                                const seissol::geometry::PumlMesh& meshGeometry) {
  bool continueComputation = true;
  if (!model::MaterialT::SupportsLTS) {
    logInfo() << "The material" << model::MaterialT::Text
              << "does not support LTS. Switching to GTS.";
    continueComputation = false;
  }
  if (rate_.empty() || (rate_.size() == 1 && rate_[0] == 1)) {
    logInfo() << "GTS has been selected.";
    continueComputation = false;
  }

  logInfo() << "Computing LTS weights.";

  // Note: Return value optimization is guaranteed while returning temp. objects in C++17
  meshTopology_ = &meshTopology;
  meshGeometry_ = &meshGeometry;
  details_ = collectGlobalTimeStepDetails();
  cellCosts_ = computeCostsPerTimestep();

  const auto& ltsParameters = seissolInstance_.getSeisSolParameters().timeStepping.lts;
  auto maxClusterIdToEnforce = ltsParameters.getMaxNumberOfClusters() - 1;

  if (!continueComputation) {
    // enforce GTS
    maxClusterIdToEnforce = 0;
  }

  prepareDifferenceEnforcement();

  if ((ltsParameters.isWiggleFactorUsed() || ltsParameters.isAutoMergeUsed()) &&
      continueComputation) {
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

    wiggleFactor_ = wiggleFactorResult.wiggleFactor;
    if (ltsParameters.isAutoMergeUsed()) {
      maxClusterIdToEnforce = std::min(maxClusterIdToEnforce, wiggleFactorResult.maxClusterId);
    }
  } else {
    wiggleFactor_ = 1.0;
  }

  ncon_ = evaluateNumberOfConstraints();
  auto finalNumberOfReductions = computeClusterIdsAndEnforceMaximumDifferenceCached(wiggleFactor_);

  if (!ltsParameters.getWiggleFactorEnforceMaximumDifference()) {
    finalNumberOfReductions += enforceMaximumDifference();
  }

  logInfo() << "Limiting number of clusters to" << maxClusterIdToEnforce + 1;
  clusterIds_ = enforceMaxClusterId(clusterIds_, maxClusterIdToEnforce);

  if (!vertexWeights_.empty()) {
    vertexWeights_.clear();
  }
  vertexWeights_.resize(clusterIds_.size() * ncon_);

  // calling virtual functions
  setVertexWeights();
  setAllowedImbalances();

  logInfo() << "Computing LTS weights. Done. " << utils::nospace << '(' << finalNumberOfReductions
            << " reductions)";
}

double LtsWeights::getWiggleFactor() const { return wiggleFactor_; }

LtsWeights::ComputeWiggleFactorResult
    LtsWeights::computeBestWiggleFactor(std::optional<double> baselineCost, bool isAutoMergeUsed) {

  // Maps that keep track of number of clusters vs cost
  auto mapMaxClusterIdToLowestCost = std::map<int, double>{};
  auto maxMapClusterIdToBestWiggleFactor = std::map<int, double>{};

  const auto& ltsParameters = seissolInstance_.getSeisSolParameters().timeStepping.lts;
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
    baselineCost = computeGlobalCostOfClustering(clusterIds_,
                                                 cellCosts_,
                                                 rate_,
                                                 maxWiggleFactor,
                                                 details_.globalMinTimeStep,
                                                 Mpi::mpi.comm());
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
          computeMaxClusterIdAfterAutoMerge(clusterIds_,
                                            cellCosts_,
                                            rate_,
                                            maxAdmissibleCost,
                                            curWiggleFactor,
                                            details_.globalMinTimeStep);
      maxClusterIdToEnforce = std::min(maxClusterIdAfterMerging, maxClusterIdToEnforce);
    }

    clusterIds_ = enforceMaxClusterId(clusterIds_, maxClusterIdToEnforce);
    auto maxClusterId = *std::max_element(clusterIds_.begin(), clusterIds_.end());
    MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, Mpi::mpi.comm());

    ncon_ = evaluateNumberOfConstraints();

    // Compute cost
    const double cost = computeGlobalCostOfClustering(clusterIds_,
                                                      cellCosts_,
                                                      rate_,
                                                      curWiggleFactor,
                                                      details_.globalMinTimeStep,
                                                      Mpi::mpi.comm());

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
  assert(!vertexWeights_.empty() && "vertex weights are not initialized");
  return vertexWeights_.data();
}

const double* LtsWeights::imbalances() const {
  assert(!imbalances_.empty() && "weight imbalances are not initialized");
  return imbalances_.data();
}

const std::vector<int>& LtsWeights::clusterIds() const { return clusterIds_; }

const std::vector<double>& LtsWeights::timesteps() const { return details_.cellTimeStepWidths; }

int LtsWeights::nWeightsPerVertex() const {
  assert(ncon_ != std::numeric_limits<int>::infinity() &&
         "num. constrains has not been initialized yet");
  return ncon_;
}

std::uint64_t LtsWeights::getCluster(double timestep,
                                     double globalMinTimestep,
                                     double ltsWiggleFactor,
                                     const std::vector<uint64_t>& rate) {
  if (rate.empty()) {
    return 0;
  }

  double upper = ltsWiggleFactor * rate[0] * globalMinTimestep;

  std::uint64_t cluster = 0;
  while (upper <= timestep) {
    const auto currentRate = rate.size() > (cluster + 1) ? rate[cluster + 1] : rate.back();
    if (currentRate == 1) {
      break;
    }
    upper *= currentRate;
    ++cluster;
  }
  return cluster;
}

FaceType LtsWeights::getBoundaryCondition(const void* boundaryCond, size_t cell, unsigned face) {
  int bcCurrentFace = seissol::geometry::decodeBoundary(boundaryCond, cell, face, boundaryFormat_);
  if (bcCurrentFace > 64) {
    bcCurrentFace = 3;
  }
  return static_cast<FaceType>(bcCurrentFace);
}

std::uint64_t ratepow(const std::vector<std::uint64_t>& rate, std::uint64_t a, std::uint64_t b) {
  std::uint64_t factor = 1;
  for (std::uint64_t i = a; i < b; ++i) {
    factor *= i < rate.size() ? rate[i] : rate.back();
  }
  return factor;
}

seissol::initializer::GlobalTimestep LtsWeights::collectGlobalTimeStepDetails() {
  return seissol::initializer::computeTimesteps(
      seissol::initializer::CellToVertexArray::fromPUML(*meshGeometry_),
      seissolInstance_.getSeisSolParameters());
}

int LtsWeights::computeClusterIdsAndEnforceMaximumDifferenceCached(double curWiggleFactor) {
  int numberOfReductions = 0;
  auto lb = clusteringCache_.lower_bound(curWiggleFactor);

  if (lb != clusteringCache_.end() && !(clusteringCache_.key_comp()(curWiggleFactor, lb->first))) {
    clusterIds_ = lb->second;
  } else {
    // re-use best computed maxdiff enforcement available
    // reason that works: cf. Lukas' proof for cluster merging not violating maximum difference
    // we may generalize due to the fact that min(a, min(b,c)) = min(min(a,b), c) = min(min(a,c),
    // b), essentially establishing a partial ordering of clusterings, where A >= B iff
    // cluster(A[i]) >= cluster(B[i]) for all cells i. Thus: walking through the wiggle factors from
    // lower to higher will save a lot of reductions

    int cellchanges = 0;
    if (lb != clusteringCache_.end()) {
      // use the cache
      const auto newClusterIds = computeClusterIds(curWiggleFactor);
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : cellchanges)
#endif
      for (std::size_t cell = 0; cell < meshTopology_->cells().size(); ++cell) {
        if (lb->second[cell] > newClusterIds[cell]) {
          ++cellchanges;
        }
        clusterIds_[cell] = std::min(lb->second[cell], newClusterIds[cell]);
      }
    } else {
      clusterIds_ = computeClusterIds(curWiggleFactor);
      cellchanges = meshTopology_->cells().size();
    }
    const auto& ltsParameters = seissolInstance_.getSeisSolParameters().timeStepping.lts;
    if (ltsParameters.getWiggleFactorEnforceMaximumDifference()) {
      MPI_Allreduce(MPI_IN_PLACE, &cellchanges, 1, MPI_INT, MPI_SUM, seissol::Mpi::mpi.comm());
      if (cellchanges > 0) {
        numberOfReductions = enforceMaximumDifference();
      }
    }
    clusteringCache_[curWiggleFactor] = clusterIds_;
  }

  return numberOfReductions;
}

std::vector<int> LtsWeights::computeClusterIds(double curWiggleFactor) {
  const auto& cells = meshTopology_->cells();
  std::vector<int> clusterIds(cells.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    clusterIds[cell] = getCluster(
        details_.cellTimeStepWidths[cell], details_.globalMinTimeStep, curWiggleFactor, rate_);
  }
  return clusterIds;
}

std::vector<int> LtsWeights::computeCostsPerTimestep() {
  const auto& cells = meshTopology_->cells();

  std::vector<int> cellCosts(cells.size());
  const void* boundaryCond = meshTopology_->cellData(1);
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    int dynamicRupture = 0;
    int freeSurfaceWithGravity = 0;

    unsigned int faceids[Cell::NumFaces];
    PUML::Downward::faces(*meshTopology_, cells[cell], faceids);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto faceType = getBoundaryCondition(boundaryCond, cell, face);
      dynamicRupture += (faceType == FaceType::DynamicRupture) ? 1 : 0;
      freeSurfaceWithGravity += (faceType == FaceType::FreeSurfaceGravity) ? 1 : 0;
    }

    const int costDynamicRupture = vertexWeightDynamicRupture_ * dynamicRupture;
    const int costDisplacement = vertexWeightFreeSurfaceWithGravity_ * freeSurfaceWithGravity;
    cellCosts[cell] = vertexWeightElement_ + costDynamicRupture + costDisplacement;
  }
  return cellCosts;
}

int LtsWeights::enforceMaximumDifference() {
  int totalNumberOfReductions = 0;
  int globalNumberOfReductions = 0;
  do {
    int localNumberOfReductions = enforceMaximumDifferenceLocal();

    MPI_Allreduce(&localNumberOfReductions,
                  &globalNumberOfReductions,
                  1,
                  MPI_INT,
                  MPI_SUM,
                  seissol::Mpi::mpi.comm());
    totalNumberOfReductions += globalNumberOfReductions;
  } while (globalNumberOfReductions > 0);
  return totalNumberOfReductions;
}

void LtsWeights::prepareDifferenceEnforcement() {
  const auto& cells = meshTopology_->cells();
  const auto& faces = meshTopology_->faces();
  const void* boundaryCond = meshTopology_->cellData(1);

  std::unordered_map<int, std::vector<std::size_t>> rankToSharedFacesPre;
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    unsigned int faceids[Cell::NumFaces]{};
    bool atBoundary = false;
    PUML::Downward::faces(*meshTopology_, cells[cell], faceids);
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbors later
        const auto& face = faces.at(faceids[f]);
        if (face.isShared()) {
          rankToSharedFacesPre[face.shared()[0]].push_back(faceids[f]);
          localFaceIdToLocalCellId_[faceids[f]] = cell;
          atBoundary = true;
        }
      }
    }
    if (atBoundary) {
      boundaryCells_.emplace_back(cell);
    }
  }

  for (auto& sharedFaces : rankToSharedFacesPre) {
    std::sort(sharedFaces.second.begin(),
              sharedFaces.second.end(),
              [&](unsigned int a, unsigned int b) { return faces[a].gid() < faces[b].gid(); });
  }

  rankToSharedFaces_ =
      decltype(rankToSharedFaces_)(rankToSharedFacesPre.begin(), rankToSharedFacesPre.end());

  for (std::size_t ex = 0; ex < rankToSharedFaces_.size(); ++ex) {
    const auto& exchange = rankToSharedFaces_[ex];
    for (std::size_t i = 0; i < exchange.second.size(); ++i) {
      sharedFaceToExchangeId_[exchange.second[i]] = {ex, i};
    }
  }
}

int LtsWeights::enforceMaximumDifferenceLocal(int maxDifference) {
  int numberOfReductions = 0;

  const auto& cells = meshTopology_->cells();
  const auto& faces = meshTopology_->faces();
  const void* boundaryCond = meshTopology_->cellData(1);

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numberOfReductions)
#endif
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    int timeCluster = clusterIds_[cell];

    unsigned int faceids[Cell::NumFaces]{};
    PUML::Downward::faces(*meshTopology_, cells[cell], faceids);
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      int difference = maxDifference;
      const auto boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbors later
        const auto& face = faces.at(faceids[f]);
        if (!face.isShared()) {
          int cellIds[2];
          PUML::Upward::cells(*meshTopology_, face, cellIds);

          const int neighborCell = (cellIds[0] == static_cast<int>(cell)) ? cellIds[1] : cellIds[0];
          const int otherTimeCluster = clusterIds_[neighborCell];

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
    clusterIds_[cell] = timeCluster;
  }

  const auto numExchanges = rankToSharedFaces_.size();
  std::vector<MPI_Request> requests(2 * numExchanges);
  std::vector<std::vector<int>> ghost(numExchanges);
  std::vector<std::vector<int>> copy(numExchanges);

  for (std::size_t ex = 0; ex < numExchanges; ++ex) {
    const auto& exchange = rankToSharedFaces_[ex];
    const auto exchangeSize = exchange.second.size();
    ghost[ex].resize(exchangeSize);
    copy[ex].resize(exchangeSize);

    for (std::size_t n = 0; n < exchangeSize; ++n) {
      copy[ex][n] = clusterIds_[localFaceIdToLocalCellId_[exchange.second[n]]];
    }
    MPI_Isend(copy[ex].data(),
              exchangeSize,
              MPI_INT,
              exchange.first,
              0,
              seissol::Mpi::mpi.comm(),
              &requests[ex]);
    MPI_Irecv(ghost[ex].data(),
              exchangeSize,
              MPI_INT,
              exchange.first,
              0,
              seissol::Mpi::mpi.comm(),
              &requests[numExchanges + ex]);
  }

  MPI_Waitall(2 * numExchanges, requests.data(), MPI_STATUSES_IGNORE);

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numberOfReductions)
#endif
  for (std::size_t bcell = 0; bcell < boundaryCells_.size(); ++bcell) {
    const auto cell = boundaryCells_[bcell];
    int& timeCluster = clusterIds_[cell];

    unsigned int faceids[Cell::NumFaces]{};
    PUML::Downward::faces(*meshTopology_, cells[cell], faceids);
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      int difference = maxDifference;
      const auto boundary = getBoundaryCondition(boundaryCond, cell, f);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbors later
        const auto& face = faces.at(faceids[f]);
        if (face.isShared()) {
          const auto pos = sharedFaceToExchangeId_.at(faceids[f]);
          const int otherTimeCluster = ghost[pos.first][pos.second];

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
  }

  return numberOfReductions;
}
} // namespace seissol::initializer::time_stepping
