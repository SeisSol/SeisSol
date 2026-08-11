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
#include "Initializer/Clustering/ClusterLadder.h"
#include "Initializer/Clustering/ClusterSmoother.h"
#include "Initializer/Clustering/GridLadderSearch.h"
#include "Initializer/Clustering/LatticeDpSearch.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"

#include <PUML/Downward.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mpi.h>
#include <optional>
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

std::vector<double> computeLocalCostsOfCappedClusterings(const std::vector<int>& clusterIds,
                                                         const std::vector<int>& cellCosts,
                                                         const std::vector<uint64_t>& rate,
                                                         double wiggleFactor,
                                                         double minimalTimestep,
                                                         int maxClusterId) {
  assert(clusterIds.size() == cellCosts.size());
  assert(maxClusterId >= 0);

  const auto capCount = static_cast<std::size_t>(maxClusterId) + 1;

  // ratepow() rather than ClusterLadder: this is also reachable with hand-made clusterings
  // whose ids exceed what the rate vector could produce, which a ladder would reject.
  std::vector<double> updateFactors(capCount);
  for (std::size_t cluster = 0; cluster < capCount; ++cluster) {
    updateFactors[cluster] = static_cast<double>(ratepow(rate, 0, cluster));
  }

  std::vector<double> costs(capCount, 0.0);
  for (auto i = 0U; i < clusterIds.size(); ++i) {
    const auto cellCost = cellCosts[i];
    const auto cluster = clusterIds[i];
    for (std::size_t cap = 0; cap < capCount; ++cap) {
      const auto capped = std::min(static_cast<std::size_t>(cluster), cap);
      costs[cap] += cellCost / updateFactors[capped];
    }
  }

  const auto minDtWithWiggle = minimalTimestep * wiggleFactor;
  for (auto& cost : costs) {
    cost /= minDtWithWiggle;
  }
  return costs;
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
  int maxClusterId =
      clusterIds.empty() ? 0 : *std::max_element(clusterIds.begin(), clusterIds.end());
  MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, Mpi::mpi.comm());

  // All candidate costs in one pass, then one reduction -- previously this was a full
  // clustering copy plus an Allreduce for every candidate.
  auto costs = computeLocalCostsOfCappedClusterings(
      clusterIds, cellCosts, rate, wiggleFactor, minimalTimestep, maxClusterId);
  MPI_Allreduce(MPI_IN_PLACE,
                costs.data(),
                static_cast<int>(costs.size()),
                MPI_DOUBLE,
                MPI_SUM,
                Mpi::mpi.comm());

  // Iteratively merge clusters until we found the first number of clusters that has a cost that is
  // too high
  for (auto curMaxClusterId = maxClusterId; curMaxClusterId >= 0; --curMaxClusterId) {
    if (costs[curMaxClusterId] > maximalAdmissibleCost) {
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
      boundaryFormat_(config.boundaryFormat), faceMap_(config.faceMap) {}

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

  evaluator_.emplace(*meshTopology_,
                     boundaryFormat_,
                     *faceMap_,
                     details_,
                     cellCosts_,
                     rate_,
                     ltsParameters.getWiggleFactorEnforceMaximumDifference());

  int finalNumberOfReductions = 0;
  std::optional<std::vector<std::uint64_t>> searchRatios;

  const auto usesLatticeSearch = ltsParameters.getClusteringSearch() ==
                                 seissol::initializer::parameters::LtsClusteringSearch::Lattice;

  if (usesLatticeSearch && ltsParameters.isAutoMergeUsed()) {
    logWarning() << "The lattice clustering search chooses the number of clusters itself;"
                 << "ltsAutoMergeClusters and its baseline setting have no effect.";
  }

  if ((ltsParameters.isWiggleFactorUsed() || ltsParameters.isAutoMergeUsed() ||
       usesLatticeSearch) &&
      continueComputation) {
    auto autoMergeBaseline = ltsParameters.getAutoMergeCostBaseline();
    if (!(ltsParameters.isWiggleFactorUsed() && ltsParameters.isAutoMergeUsed())) {
      // Cost models only change things if both wiggle factor and auto merge are on.
      // In all other cases, choose the cheapest cost model.
      autoMergeBaseline = seissol::initializer::parameters::AutoMergeCostBaseline::MaxWiggleFactor;
    }

    const SearchConstraints constraints{ltsParameters.getWiggleFactorMinimum(),
                                        ltsParameters.getWiggleFactorStepsize(),
                                        ltsParameters.getMaxNumberOfClusters() - 1,
                                        ltsParameters.isAutoMergeUsed(),
                                        ltsParameters.getAllowedPerformanceLossRatioAutoMerge(),
                                        autoMergeBaseline,
                                        ltsParameters.getClusterCostModel(),
                                        ltsParameters.getMaxRate()};

    std::unique_ptr<LadderSearch> search;
    if (usesLatticeSearch) {
      search = std::make_unique<LatticeDpSearch>();
    } else {
      search = std::make_unique<GridLadderSearch>();
    }
    const auto searchResult = search->run(evaluator_.value(), constraints);

    wiggleFactor_ = searchResult.wiggleFactor;
    searchRatios = searchResult.ratios;
    if (ltsParameters.isAutoMergeUsed()) {
      maxClusterIdToEnforce = std::min(maxClusterIdToEnforce, searchResult.maxClusterId);
    }
  } else {
    wiggleFactor_ = 1.0;
  }

  // Normalize once, here: from this point on nothing may re-derive the ladder from the
  // parameter file, because the search is allowed to have changed it.
  effectiveRates_ = searchRatios.value_or(
      ClusterLadder::forBinning(
          rate_, details_.globalMinTimeStep, wiggleFactor_, details_.globalMaxTimeStep)
          .ratios());

  ncon_ = evaluateNumberOfConstraints();
  finalNumberOfReductions += evaluator_.value().realize(effectiveRates_, wiggleFactor_);

  if (!ltsParameters.getWiggleFactorEnforceMaximumDifference()) {
    finalNumberOfReductions += evaluator_.value().smoothCurrent();
  }

  clusterIds_ = evaluator_.value().clusterIds();

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
  logInfo() << "Cluster rates:" << effectiveRates_;
}

double LtsWeights::getWiggleFactor() const { return wiggleFactor_; }

const std::vector<std::uint64_t>& LtsWeights::effectiveRates() const { return effectiveRates_; }

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

std::uint64_t getCluster(double timestep,
                         double globalMinTimestep,
                         double ltsWiggleFactor,
                         const std::vector<std::uint64_t>& rate) {
  // Kept as a standalone loop rather than going through ClusterLadder: this is called for
  // single ad-hoc queries where building a ladder would allocate. It shares the ladder's
  // notion of where the rate vector terminates, and the unit tests pin the two together.
  if (rate.empty()) {
    return 0;
  }

  const auto clusterLimit = ClusterLadder::intrinsicClusterCount(rate);

  double upper = ltsWiggleFactor * rate[0] * globalMinTimestep;

  std::uint64_t cluster = 0;
  while (upper <= timestep) {
    if (clusterLimit != ClusterLadder::Unbounded && cluster + 1 >= clusterLimit) {
      break;
    }
    const auto currentRate = rate.size() > (cluster + 1) ? rate[cluster + 1] : rate.back();
    upper *= currentRate;
    ++cluster;
  }
  return cluster;
}

FaceType LtsWeights::getBoundaryCondition(const void* boundaryCond, size_t cell, unsigned face) {
  return decodeFaceType(boundaryCond, cell, face, boundaryFormat_, *faceMap_);
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

} // namespace seissol::initializer::time_stepping
