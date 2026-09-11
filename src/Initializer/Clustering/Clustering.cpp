// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/Clustering/Clustering.h"

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Clustering/ClusterLadder.h"
#include "Initializer/Clustering/ClusterSmoother.h"
#include "Initializer/Clustering/ClusteringCost.h"
#include "Initializer/Clustering/ClusteringEvaluator.h"
#include "Initializer/Clustering/GridLadderSearch.h"
#include "Initializer/Clustering/LadderSearch.h"
#include "Initializer/Clustering/LatticeDpSearch.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/LtsParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "SeisSol.h"

#include <PUML/Downward.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

Clustering::Clustering(const ClusteringConfig& config, seissol::SeisSol& seissolInstance)
    : seissolInstance_(seissolInstance), rate_(config.rate),
      vertexWeightElement_(config.vertexWeightElement),
      vertexWeightDynamicRupture_(config.vertexWeightDynamicRupture),
      vertexWeightFreeSurfaceWithGravity_(config.vertexWeightFreeSurfaceWithGravity),
      boundaryFormat_(config.boundaryFormat), faceMap_(config.faceMap) {}

const ClusteringResult& Clustering::compute(const seissol::geometry::PumlMesh& meshTopology,
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

  auto details =
      computeTimesteps(CellToVertexArray::fromPUML(meshGeometry), seissolInstance_.parameters());
  auto cellCosts = computeCostsPerTimestep(meshTopology);

  const auto& ltsParameters = seissolInstance_.parameters().timeStepping.lts;
  // getMaxNumberOfClusters() is validated to be positive, so the decrement cannot wrap
  auto maxClusterIdToEnforce = static_cast<std::size_t>(ltsParameters.getMaxNumberOfClusters()) - 1;

  if (!continueComputation) {
    // enforce GTS
    maxClusterIdToEnforce = 0;
  }

  // Scoped: the evaluator borrows `details` and `cellCosts`, which are moved into the result
  // once it is gone.
  ClusteringEvaluator evaluator(meshTopology,
                                boundaryFormat_,
                                *faceMap_,
                                details,
                                cellCosts,
                                rate_,
                                ltsParameters.getWiggleFactorEnforceMaximumDifference());

  std::size_t finalNumberOfReductions = 0;
  double wiggleFactor = 1.0;
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
                                        maxClusterIdToEnforce,
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
    const auto searchResult = search->run(evaluator, constraints);

    wiggleFactor = searchResult.wiggleFactor;
    searchRatios = searchResult.ratios;
    if (ltsParameters.isAutoMergeUsed()) {
      maxClusterIdToEnforce = std::min(maxClusterIdToEnforce, searchResult.maxClusterId);
    }
  }

  // Normalize once, here: from this point on nothing may re-derive the ladder from the
  // parameter file, because the search is allowed to have changed it.
  auto effectiveRates = searchRatios.value_or(
      ClusterLadder::forBinning(
          rate_, details.globalMinTimeStep, wiggleFactor, details.globalMaxTimeStep)
          .ratios());

  finalNumberOfReductions += evaluator.realize(effectiveRates, wiggleFactor);

  if (!ltsParameters.getWiggleFactorEnforceMaximumDifference()) {
    finalNumberOfReductions += evaluator.smoothCurrent();
  }

  auto clusterIds = evaluator.clusterIds();

  logInfo() << "Limiting number of clusters to" << maxClusterIdToEnforce + 1;
  clusterIds = enforceMaxClusterId(clusterIds, maxClusterIdToEnforce);

  logInfo() << "Computing LTS weights. Done. " << utils::nospace << '(' << finalNumberOfReductions
            << " reductions)";
  logInfo() << "Cluster rates:" << effectiveRates;

  result_ = ClusteringResult{std::move(clusterIds),
                             std::move(effectiveRates),
                             wiggleFactor,
                             std::move(cellCosts),
                             std::move(details)};
  return result_;
}

std::vector<std::uint64_t>
    Clustering::computeCostsPerTimestep(const seissol::geometry::PumlMesh& mesh) const {
  const auto& cells = mesh.cells();

  std::vector<std::uint64_t> cellCosts(cells.size());
  const void* boundaryCond = mesh.cellData(1);
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    std::uint64_t dynamicRupture = 0;
    std::uint64_t freeSurfaceWithGravity = 0;

    unsigned int faceids[Cell::NumFaces];
    PUML::Downward::faces(mesh, cells[cell], faceids);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto faceType = decodeFaceType(boundaryCond, cell, face, boundaryFormat_, *faceMap_);
      dynamicRupture += (faceType == FaceType::DynamicRupture) ? 1 : 0;
      freeSurfaceWithGravity += (faceType == FaceType::FreeSurfaceGravity) ? 1 : 0;
    }

    const auto costDynamicRupture = vertexWeightDynamicRupture_ * dynamicRupture;
    const auto costDisplacement = vertexWeightFreeSurfaceWithGravity_ * freeSurfaceWithGravity;
    cellCosts[cell] = vertexWeightElement_ + costDynamicRupture + costDisplacement;
  }
  return cellCosts;
}

} // namespace seissol::initializer
