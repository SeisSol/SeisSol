// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "Initializer/Clustering/ClusteringEvaluator.h"

#include "Geometry/PUMLReader.h"
#include "Initializer/Clustering/ClusterHistogram.h"
#include "Initializer/Clustering/ClusterLadder.h"
#include "Initializer/Clustering/ClusterSmoother.h"
#include "Initializer/Clustering/TimestepHistogram.h"
#include "Initializer/Clustering/VertexWeights/LtsWeights.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <utility>
#include <vector>

namespace {
using seissol::initializer::time_stepping::computeGlobalCostOfClustering;
} // namespace

namespace seissol::initializer {

ClusteringEvaluator::ClusteringEvaluator(const geometry::PumlMesh& mesh,
                                         parameters::BoundaryFormat boundaryFormat,
                                         const FaceMap& faceMap,
                                         const GlobalTimestep& timesteps,
                                         const std::vector<int>& cellCosts,
                                         std::vector<std::uint64_t> rate,
                                         bool smoothDuringSearch)
    : smoother_(mesh, boundaryFormat, faceMap), timesteps_(&timesteps), cellCosts_(&cellCosts),
      rate_(std::move(rate)), smoothDuringSearch_(smoothDuringSearch),
      cellCount_(mesh.cells().size()), clusterIds_(mesh.cells().size(), 0) {}

ClusterLadder ClusteringEvaluator::configuredLadder(double wiggleFactor) const {
  // The only place that applies the abbreviated ClusteredLTS convention; everything past the
  // search works on complete ladders.
  return ClusterLadder::forBinning(
      rate_, timesteps_->globalMinTimeStep, wiggleFactor, timesteps_->globalMaxTimeStep);
}

std::vector<std::uint64_t> ClusteringEvaluator::configuredRatios(double wiggleFactor) const {
  return configuredLadder(wiggleFactor).ratios();
}

std::vector<int> ClusteringEvaluator::binCells(const ClusterLadder& ladder) const {
  std::vector<int> clusterIds(cellCount_, 0);

#pragma omp parallel for
  for (std::size_t cell = 0; cell < cellCount_; ++cell) {
    clusterIds[cell] = static_cast<int>(ladder.clusterOf(timesteps_->cellTimeStepWidths[cell]));
  }
  return clusterIds;
}

int ClusteringEvaluator::realize(double wiggleFactor) {
  int numberOfReductions = 0;
  auto lb = cache_.lower_bound(wiggleFactor);

  if (lb != cache_.end() && !(cache_.key_comp()(wiggleFactor, lb->first))) {
    clusterIds_ = lb->second;
  } else {
    // re-use best computed maxdiff enforcement available
    // reason that works: cf. Lukas' proof for cluster merging not violating maximum difference
    // we may generalize due to the fact that min(a, min(b,c)) = min(min(a,b), c) = min(min(a,c),
    // b), essentially establishing a partial ordering of clusterings, where A >= B iff
    // cluster(A[i]) >= cluster(B[i]) for all cells i. Thus: walking through the wiggle factors from
    // lower to higher will save a lot of reductions

    int cellchanges = 0;
    if (lb != cache_.end()) {
      // use the cache
      const auto newClusterIds = binCells(configuredLadder(wiggleFactor));

#pragma omp parallel for reduction(+ : cellchanges)
      for (std::size_t cell = 0; cell < cellCount_; ++cell) {
        if (lb->second[cell] > newClusterIds[cell]) {
          ++cellchanges;
        }
        clusterIds_[cell] = std::min(lb->second[cell], newClusterIds[cell]);
      }
    } else {
      clusterIds_ = binCells(configuredLadder(wiggleFactor));
      cellchanges = static_cast<int>(cellCount_);
    }
    if (smoothDuringSearch_) {
      MPI_Allreduce(MPI_IN_PLACE, &cellchanges, 1, MPI_INT, MPI_SUM, seissol::Mpi::mpi.comm());
      if (cellchanges > 0) {
        numberOfReductions = smoothCurrent();
      }
    }
    cache_[wiggleFactor] = clusterIds_;
  }

  return numberOfReductions;
}

int ClusteringEvaluator::realize(const std::vector<std::uint64_t>& ratios, double wiggleFactor) {
  if (ratios == configuredRatios(wiggleFactor)) {
    return realize(wiggleFactor);
  }
  clusterIds_ = binCells(ClusterLadder::exact(ratios, timesteps_->globalMinTimeStep, wiggleFactor));
  return smoothDuringSearch_ ? smoothCurrent() : 0;
}

int ClusteringEvaluator::smoothCurrent() {
  return smoother_.relax(clusterIds_, smoothingRule_, seissol::Mpi::mpi.comm());
}

int ClusteringEvaluator::globalMaxClusterId() const {
  int maxClusterId =
      clusterIds_.empty() ? 0 : *std::max_element(clusterIds_.begin(), clusterIds_.end());
  MPI_Allreduce(MPI_IN_PLACE, &maxClusterId, 1, MPI_INT, MPI_MAX, Mpi::mpi.comm());
  return maxClusterId;
}

double ClusteringEvaluator::globalCost(double wiggleFactor) const {
  return computeGlobalCostOfClustering(clusterIds_,
                                       *cellCosts_,
                                       rate_,
                                       wiggleFactor,
                                       timesteps_->globalMinTimeStep,
                                       Mpi::mpi.comm());
}

double ClusteringEvaluator::globalCost(const std::vector<std::uint64_t>& ratios,
                                       double wiggleFactor) const {
  return computeGlobalCostOfClustering(clusterIds_,
                                       *cellCosts_,
                                       ratios,
                                       wiggleFactor,
                                       timesteps_->globalMinTimeStep,
                                       Mpi::mpi.comm());
}

TimestepHistogram ClusteringEvaluator::timestepHistogram(double wiggleFactor,
                                                         std::size_t maxIndex) const {
  auto histogram = TimestepHistogram::fromCells(timesteps_->cellTimeStepWidths,
                                                *cellCosts_,
                                                timesteps_->globalMinTimeStep * wiggleFactor,
                                                maxIndex);
  histogram.reduce(Mpi::mpi.comm());
  return histogram;
}

ClusterHistogram ClusteringEvaluator::globalHistogram() const {
  auto histogram = ClusterHistogram::fromClustering(
      clusterIds_, *cellCosts_, static_cast<std::size_t>(globalMaxClusterId()) + 1);
  histogram.reduce(Mpi::mpi.comm());
  return histogram;
}

} // namespace seissol::initializer
