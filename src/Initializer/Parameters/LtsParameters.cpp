// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LtsParameters.h"

#include "Equations/Datastructures.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::initializer::parameters {

AutoMergeCostBaseline parseAutoMergeCostBaseline(std::string str) {
  // Convert str to lower case to make function case-insensitive
  // Note: This is, of course, broken for non-ASCI input.
  std::transform(str.begin(), str.end(), str.begin(), [](auto c) { return std::tolower(c); });
  if (str == "bestwigglefactor") {
    return AutoMergeCostBaseline::BestWiggleFactor;
  } else if (str == "maxwigglefactor") {
    return AutoMergeCostBaseline::MaxWiggleFactor;
  }
  throw std::invalid_argument(str + " is not a valid cluster merging baseline");
}

LtsClusteringSearch parseLtsClusteringSearch(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), [](auto c) { return std::tolower(c); });
  if (str == "legacy") {
    return LtsClusteringSearch::Legacy;
  } else if (str == "lattice") {
    return LtsClusteringSearch::Lattice;
  }
  throw std::invalid_argument(str + " is not a valid LTS clustering search");
}

LtsParameters readLtsParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("discretization");
  const auto ratestr = reader->readWithDefault<std::string>("clusteredlts", "1");
  std::vector<uint64_t> rates;
  auto parts = utils::StringUtils::split(ratestr, ' ');
  for (auto& part : parts) {
    utils::StringUtils::trim(part);
    rates.emplace_back(std::stoull(part));

    if (rates.back() == 0) {
      logError() << "Invalid LTS rate (0) found in" << ratestr << "after parsing" << rates
                 << ". Aborting.";
    }
  }

  if (rates.empty()) {
    logWarning() << "No LTS rate given. Assuming GTS.";
    rates.emplace_back(1);
  }

  for (std::size_t i = 0; i + 1 < rates.size(); ++i) {
    if (rates[i] == 1) {
      logError() << "Invalid LTS rate (1) found in" << rates << ". Aborting.";
    }
  }

  const double wiggleFactorMinimum = reader->readWithDefault("ltswigglefactormin", 1.0);
  const double wiggleFactorStepsize = reader->readWithDefault("ltswigglefactorstepsize", 0.01);
  const bool wiggleFactorEnforceMaximumDifference =
      reader->readWithDefault("ltswigglefactorenforcemaximumdifference", true);
  const unsigned int maxNumberOfClusters =
      reader->readWithDefault("ltsmaxnumberofclusters", std::numeric_limits<int>::max() - 1);
  const bool autoMergeClusters = reader->readWithDefault("ltsautomergeclusters", false);
  const double allowedRelativePerformanceLossAutoMerge =
      reader->readWithDefault("ltsallowedrelativeperformancelossautomerge", 0.0);
  const double allowedPerformanceLossRatioAutoMerge = allowedRelativePerformanceLossAutoMerge + 1.0;
  const auto autoMergeCostBaseline = parseAutoMergeCostBaseline(
      (reader->readWithDefault("ltsautomergecostbaseline", std::string("bestwigglefactor"))));
  const LtsWeightsTypes ltsWeightsType =
      reader->readWithDefaultEnum("ltsweighttypeid",
                                  LtsWeightsTypes::ExponentialWeights,
                                  {
                                      LtsWeightsTypes::ExponentialWeights,
                                      LtsWeightsTypes::ExponentialBalancedWeights,
                                      LtsWeightsTypes::EncodedBalancedWeights,
                                  });

  const auto clusteringSearch = parseLtsClusteringSearch(
      reader->readWithDefault("ltsclusteringsearch", std::string("legacy")));
  // Defaults reproduce the update-count objective SeisSol has always minimized; see
  // ClusterCostModel.
  const double costUpdate = reader->readWithDefault("ltscostupdate", 1.0);
  const double costLaunch = reader->readWithDefault("ltscostlaunch", 0.0);
  const double costFill = reader->readWithDefault("ltscostfill", 0.0);
  const auto maxRate = static_cast<std::uint64_t>(reader->readWithDefault("ltsmaxrate", 0));

  reader->warnDeprecated({"dgmethod"});
  return {rates,
          wiggleFactorMinimum,
          wiggleFactorStepsize,
          wiggleFactorEnforceMaximumDifference,
          static_cast<int>(maxNumberOfClusters),
          autoMergeClusters,
          allowedPerformanceLossRatioAutoMerge,
          autoMergeCostBaseline,
          ltsWeightsType,
          clusteringSearch,
          costUpdate,
          costLaunch,
          costFill,
          maxRate};
}

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
  return ClusterCostModel{costUpdate_, costLaunch_, costFill_};
}

TimeSteppingParameters::TimeSteppingParameters(VertexWeightParameters vertexWeight,
                                               double cfl,
                                               double maxTimestepWidth,
                                               double endTime,
                                               LtsParameters lts)
    : vertexWeight(vertexWeight), cfl(cfl), maxTimestepWidth(maxTimestepWidth), endTime(endTime),
      lts(std::move(lts)) {}

TimeSteppingParameters readTimeSteppingParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("discretization");
  const auto weightElement = reader->readWithDefault("vertexweightelement", 100);
  const auto weightDynamicRupture = reader->readWithDefault("vertexweightdynamicrupture", 100);
  const auto weightFreeSurfaceWithGravity =
      reader->readWithDefault("vertexweightfreesurfacewithgravity", 100);
  const double cfl = reader->readWithDefault("cfl", 0.5);
  double maxTimestepWidth = std::numeric_limits<double>::max();

  constexpr auto IsAnelastic = seissol::model::MaterialT::Mechanisms > 0;

  if constexpr (IsAnelastic) {
    auto* modelReader = baseReader->readSubNode("equations");
    const auto freqCentral = modelReader->readIfRequired<double>("freqcentral", IsAnelastic);
    const auto freqRatio = modelReader->readIfRequired<double>("freqratio", IsAnelastic);
    const double maxTimestepWidthDefault = 0.25 / (freqCentral * std::sqrt(freqRatio));
    maxTimestepWidth = reader->readWithDefault("fixtimestep", maxTimestepWidthDefault);
    if (maxTimestepWidth > maxTimestepWidthDefault) {
      logWarning()
          << "The given maximum timestep width (fixtimestep) is set to" << maxTimestepWidth
          << "which is larger than the recommended value of" << maxTimestepWidthDefault
          << "for visco-elastic material (as specified in the documentation). Please be aware"
             "that a too large maximum timestep width may cause the solution to become unstable.";
    } else {
      logInfo() << "Maximum timestep width (fixtimestep) given as" << maxTimestepWidth
                << "(less or equal to reference timestep" << maxTimestepWidthDefault << ")";
    }
  } else {
    maxTimestepWidth = reader->readWithDefault("fixtimestep", 5000.0);
  }

  auto* timeReader = baseReader->readSubNode("abortcriteria");
  const double endTime = timeReader->readWithDefault("endtime", 15.0);

  const LtsParameters ltsParameters = readLtsParameters(baseReader);

  reader->warnDeprecated({"ckmethod",
                          "dgfineout1d",
                          "fluxmethod",
                          "iterationcriterion",
                          "npoly",
                          "npolyrec",
                          "limitersecurityfactor",
                          "order",
                          "material",
                          "npolymap"});

  return TimeSteppingParameters({weightElement, weightDynamicRupture, weightFreeSurfaceWithGravity},
                                cfl,
                                maxTimestepWidth,
                                endTime,
                                ltsParameters);
}

} // namespace seissol::initializer::parameters
