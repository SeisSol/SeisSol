// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LtsParameters.h"

#include "Core/TimeStepping/Clustering/ClusterCostModel.h"

#include <cstdint>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::parameters {

LtsParameters::LtsParameters(const std::vector<uint64_t>& rates,
                             double wiggleFactorMinimum,
                             double wiggleFactorStepsize,
                             bool wigleFactorEnforceMaximumDifference,
                             int maxNumberOfClusters,
                             bool ltsAutoMergeClusters,
                             double allowedPerformanceLossRatioAutoMerge,
                             AutoMergeCostBaseline autoMergeCostBaseline,
                             LtsWeightsTypes ltsWeightsType,
                             LtsClusteringSearch clusteringSearch,
                             double costUpdate,
                             double costLaunch,
                             double costFill,
                             std::uint64_t maxRate)
    : rate_(rates), wiggleFactorMinimum_(wiggleFactorMinimum),
      wiggleFactorStepsize_(wiggleFactorStepsize),
      wiggleFactorEnforceMaximumDifference_(wigleFactorEnforceMaximumDifference),
      maxNumberOfClusters_(maxNumberOfClusters), autoMergeClusters_(ltsAutoMergeClusters),
      allowedPerformanceLossRatioAutoMerge_(allowedPerformanceLossRatioAutoMerge),
      autoMergeCostBaseline_(autoMergeCostBaseline), ltsWeightsType_(ltsWeightsType),
      clusteringSearch_(clusteringSearch), costUpdate_(costUpdate), costLaunch_(costLaunch),
      costFill_(costFill), maxRate_(maxRate) {

  if (rate_.empty()) {
    rate_.emplace_back(1);
  }

  const bool isWiggleFactorValid =
      (rate_[0] == 1 && wiggleFactorMinimum == 1.0) ||
      (wiggleFactorMinimum <= 1.0 && wiggleFactorMinimum > (1.0 / rate_[0]));
  if (!isWiggleFactorValid) {
    logError() << "Minimal wiggle factor of " << wiggleFactorMinimum << "is not valid for rate"
               << rate_;
  }
  if (maxNumberOfClusters <= 0) {
    logError() << "At least one cluster is required. Settings ltsMaxNumberOfClusters is invalid.";
  }
  if (allowedPerformanceLossRatioAutoMerge < 1.0) {
    logError() << "Negative performance loss for auto merge is invalid.";
  }
  if (costUpdate_ <= 0.0) {
    logError() << "The LTS update cost has to be positive.";
  }
  if (costLaunch_ < 0.0 || costFill_ < 0.0) {
    logError() << "The LTS launch cost and fill threshold must not be negative.";
  }
  if (maxRate_ == 1) {
    logError() << "A maximum LTS rate of 1 would forbid any clustering.";
  }
  if (clusteringSearch_ != LtsClusteringSearch::Lattice && maxRate_ != 0) {
    logWarning() << "ltsMaxRate only has an effect for the lattice clustering search; ignoring it.";
  }
}

bool LtsParameters::isWiggleFactorUsed() const { return wiggleFactorMinimum_ < 1.0; }

std::vector<uint64_t> LtsParameters::getRate() const { return rate_; }

LtsWeightsTypes LtsParameters::getLtsWeightsType() const { return ltsWeightsType_; }

double LtsParameters::getWiggleFactorMinimum() const { return wiggleFactorMinimum_; }

double LtsParameters::getWiggleFactorStepsize() const { return wiggleFactorStepsize_; }

bool LtsParameters::getWiggleFactorEnforceMaximumDifference() const {
  return wiggleFactorEnforceMaximumDifference_;
}

int LtsParameters::getMaxNumberOfClusters() const { return maxNumberOfClusters_; }

bool LtsParameters::isAutoMergeUsed() const { return autoMergeClusters_; }

double LtsParameters::getAllowedPerformanceLossRatioAutoMerge() const {
  return allowedPerformanceLossRatioAutoMerge_;
}
AutoMergeCostBaseline LtsParameters::getAutoMergeCostBaseline() const {
  return autoMergeCostBaseline_;
}

LtsClusteringSearch LtsParameters::getClusteringSearch() const { return clusteringSearch_; }

std::uint64_t LtsParameters::getMaxRate() const { return maxRate_; }

ClusterCostModel LtsParameters::getClusterCostModel() const {
  return ClusterCostModel{
      .updateCost = costUpdate_, .launchCost = costLaunch_, .fillThreshold = costFill_};
}

TimeSteppingParameters::TimeSteppingParameters(VertexWeightParameters vertexWeight,
                                               double cfl,
                                               double maxTimestepWidth,
                                               double endTime,
                                               LtsParameters lts)
    : vertexWeight(std::move(vertexWeight)), cfl(cfl), maxTimestepWidth(maxTimestepWidth),
      endTime(endTime), lts(std::move(lts)) {}

} // namespace seissol::initializer::parameters
