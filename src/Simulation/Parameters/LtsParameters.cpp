// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LtsParameters.h"

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/ConfigValue.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "Simulation/Parameters/ParameterReader.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::initializer::parameters {

AutoMergeCostBaseline parseAutoMergeCostBaseline(std::string str) {
  // Convert str to lower case to make function case-insensitive
  // Note: This is, of course, broken for non-ASCI input.
  std::ranges::transform(str, str.begin(), [](auto c) { return std::tolower(c); });
  if (str == "bestwigglefactor") {
    return AutoMergeCostBaseline::BestWiggleFactor;
  } else if (str == "maxwigglefactor") {
    return AutoMergeCostBaseline::MaxWiggleFactor;
  }
  throw std::invalid_argument(str + " is not a valid cluster merging baseline");
}

LtsClusteringSearch parseLtsClusteringSearch(std::string str) {
  std::ranges::transform(str, str.begin(), [](auto c) { return std::tolower(c); });
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

TimeSteppingParameters readTimeSteppingParameters(ParameterReader* baseReader,
                                                  const std::vector<ConfigId>& configs) {
  auto* reader = baseReader->readSubNode("discretization");
  const auto weightElement =
      static_cast<std::uint64_t>(reader->readWithDefault("vertexweightelement", 100));
  const auto weightDynamicRupture =
      static_cast<std::uint64_t>(reader->readWithDefault("vertexweightdynamicrupture", 100));
  const auto weightFreeSurfaceWithGravity = static_cast<std::uint64_t>(
      reader->readWithDefault("vertexweightfreesurfacewithgravity", 100));
  const auto measuredCosts = reader->readWithDefault<std::string>("vertexweightfile", "");
  const auto cycles = reader->readWithDefault<unsigned>("vertexweightcycles", 0);
  const auto cycleSteps = reader->readWithDefault<unsigned>("vertexweightcyclesteps", 20);
  if (cycles > 0 && cycleSteps == 0) {
    logError() << "vertexweightcyclesteps has to be positive.";
  }
  const double cfl = reader->readWithDefault("cfl", 0.5);
  double maxTimestepWidth = std::numeric_limits<double>::max();

  const bool isAnelastic = std::ranges::any_of(
      configs, [](ConfigId config) { return configValue(config).relaxationMechanisms > 0; });

  if (isAnelastic) {
    auto* modelReader = baseReader->readSubNode("equations");
    const auto freqCentral = modelReader->readIfRequired<double>("freqcentral", isAnelastic);
    const auto freqRatio = modelReader->readIfRequired<double>("freqratio", isAnelastic);
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

  return TimeSteppingParameters({.weightElement = weightElement,
                                 .weightDynamicRupture = weightDynamicRupture,
                                 .weightFreeSurfaceWithGravity = weightFreeSurfaceWithGravity,
                                 .measuredCosts = measuredCosts,
                                 .cycles = cycles,
                                 .cycleSteps = cycleSteps},
                                cfl,
                                maxTimestepWidth,
                                endTime,
                                ltsParameters);
}

} // namespace seissol::initializer::parameters
