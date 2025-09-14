// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LtsParameters.h"

#include <Equations/Datastructures.h>
#include <Initializer/Parameters/ParameterReader.h>
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
  return {rates,
          wiggleFactorMinimum,
          wiggleFactorStepsize,
          wiggleFactorEnforceMaximumDifference,
          static_cast<int>(maxNumberOfClusters),
          autoMergeClusters,
          allowedPerformanceLossRatioAutoMerge,
          autoMergeCostBaseline,
          ltsWeightsType};
}

LtsParameters::LtsParameters(const std::vector<uint64_t>& rates,
                             double wiggleFactorMinimum,
                             double wiggleFactorStepsize,
                             bool wigleFactorEnforceMaximumDifference,
                             int maxNumberOfClusters,
                             bool ltsAutoMergeClusters,
                             double allowedPerformanceLossRatioAutoMerge,
                             AutoMergeCostBaseline autoMergeCostBaseline,
                             LtsWeightsTypes ltsWeightsType)
    : rate(rates), wiggleFactorMinimum(wiggleFactorMinimum),
      wiggleFactorStepsize(wiggleFactorStepsize),
      wiggleFactorEnforceMaximumDifference(wigleFactorEnforceMaximumDifference),
      maxNumberOfClusters(maxNumberOfClusters), autoMergeClusters(ltsAutoMergeClusters),
      allowedPerformanceLossRatioAutoMerge(allowedPerformanceLossRatioAutoMerge),
      autoMergeCostBaseline(autoMergeCostBaseline), ltsWeightsType(ltsWeightsType) {

  if (rate.empty()) {
    rate.emplace_back(1);
  }

  const bool isWiggleFactorValid =
      (rate[0] == 1 && wiggleFactorMinimum == 1.0) ||
      (wiggleFactorMinimum <= 1.0 && wiggleFactorMinimum > (1.0 / rate[0]));
  if (!isWiggleFactorValid) {
    logError() << "Minimal wiggle factor of " << wiggleFactorMinimum << "is not valid for rate"
               << rate;
  }
  if (maxNumberOfClusters <= 0) {
    logError() << "At least one cluster is required. Settings ltsMaxNumberOfClusters is invalid.";
  }
  if (allowedPerformanceLossRatioAutoMerge < 1.0) {
    logError() << "Negative performance loss for auto merge is invalid.";
  }
}

bool LtsParameters::isWiggleFactorUsed() const { return wiggleFactorMinimum < 1.0; }

std::vector<uint64_t> LtsParameters::getRate() const { return rate; }

LtsWeightsTypes LtsParameters::getLtsWeightsType() const { return ltsWeightsType; }

double LtsParameters::getWiggleFactorMinimum() const { return wiggleFactorMinimum; }

double LtsParameters::getWiggleFactorStepsize() const { return wiggleFactorStepsize; }

double LtsParameters::getWiggleFactor() const { return finalWiggleFactor; }

bool LtsParameters::getWiggleFactorEnforceMaximumDifference() const {
  return wiggleFactorEnforceMaximumDifference;
}

int LtsParameters::getMaxNumberOfClusters() const { return maxNumberOfClusters; }

bool LtsParameters::isAutoMergeUsed() const { return autoMergeClusters; }

double LtsParameters::getAllowedPerformanceLossRatioAutoMerge() const {
  return allowedPerformanceLossRatioAutoMerge;
}
AutoMergeCostBaseline LtsParameters::getAutoMergeCostBaseline() const {
  return autoMergeCostBaseline;
}

void LtsParameters::setWiggleFactor(double factor) {
  assert(factor >= 1.0 / static_cast<double>(rate[0]));
  assert(factor <= 1.0);
  finalWiggleFactor = factor;
}

void LtsParameters::setMaxNumberOfClusters(int numClusters) {
  assert(numClusters > 0);
  maxNumberOfClusters = numClusters;
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
