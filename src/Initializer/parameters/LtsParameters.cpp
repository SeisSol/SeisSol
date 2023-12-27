#include "LtsParameters.h"
#include <utils/logger.h>

namespace seissol::initializers::parameters {

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

LtsParameters readLtsParameters(ParameterReader& baseReader) {
  auto reader = baseReader.readSubNode("discretization");
  const unsigned int rate = reader.readWithDefault("clusteredlts", 1);
  const double wiggleFactorMinimum = reader.readWithDefault("ltswigglefactormin", 1.0);
  const double wiggleFactorStepsize = reader.readWithDefault("ltswigglefactorstepsize", 0.01);
  const bool wiggleFactorEnforceMaximumDifference =
      reader.readWithDefault("ltswigglefactorenforcemaximumdifference", true);
  const unsigned int maxNumberOfClusters =
      reader.readWithDefault("ltsmaxnumberofclusters", std::numeric_limits<int>::max() - 1);
  const bool autoMergeClusters = reader.readWithDefault("ltsautomergeclusters", false);
  const double allowedRelativePerformanceLossAutoMerge =
      reader.readWithDefault("ltsallowedrelativeperformancelossautomerge", 0.0);
  const double allowedPerformanceLossRatioAutoMerge = allowedRelativePerformanceLossAutoMerge + 1.0;
  const auto autoMergeCostBaseline = parseAutoMergeCostBaseline(
      (reader.readWithDefault("ltsautomergecostbaseline", std::string("bestwigglefactor"))));
  const LtsWeightsTypes ltsWeightsType =
      reader.readWithDefaultEnum("ltsweighttypeid",
                                 LtsWeightsTypes::ExponentialWeights,
                                 {
                                     LtsWeightsTypes::ExponentialWeights,
                                     LtsWeightsTypes::ExponentialBalancedWeights,
                                     LtsWeightsTypes::EncodedBalancedWeights,
                                 });
  return LtsParameters(rate,
                       wiggleFactorMinimum,
                       wiggleFactorStepsize,
                       wiggleFactorEnforceMaximumDifference,
                       maxNumberOfClusters,
                       autoMergeClusters,
                       allowedPerformanceLossRatioAutoMerge,
                       autoMergeCostBaseline);
}

LtsParameters::LtsParameters(unsigned int rate,
                             double wiggleFactorMinimum,
                             double wiggleFactorStepsize,
                             bool wigleFactorEnforceMaximumDifference,
                             int maxNumberOfClusters,
                             bool ltsAutoMergeClusters,
                             double allowedPerformanceLossRatioAutoMerge,
                             AutoMergeCostBaseline autoMergeCostBaseline)
    : rate(rate), wiggleFactorMinimum(wiggleFactorMinimum),
      wiggleFactorStepsize(wiggleFactorStepsize),
      wiggleFactorEnforceMaximumDifference(wigleFactorEnforceMaximumDifference),
      maxNumberOfClusters(maxNumberOfClusters), autoMergeClusters(ltsAutoMergeClusters),
      allowedPerformanceLossRatioAutoMerge(allowedPerformanceLossRatioAutoMerge),
      autoMergeCostBaseline(autoMergeCostBaseline) {
  const bool isWiggleFactorValid =
      (rate == 1 && wiggleFactorMinimum == 1.0) ||
      (wiggleFactorMinimum <= 1.0 && wiggleFactorMinimum > (1.0 / rate));
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

unsigned int LtsParameters::getRate() const { return rate; }

LtsWeightsTypes LtsParameters::getLtsWeightsType() const { return ltsWeightsType; }

double LtsParameters::getWiggleFactorMinimum() const { return wiggleFactorMinimum; }

double LtsParameters::getWiggleFactorStepsize() const { return wiggleFactorStepsize; }

double LtsParameters::getWiggleFactor() const { return wiggleFactor; }

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

TimeSteppingParameters::TimeSteppingParameters(VertexWeightParameters vertexWeight,
                                               double cfl,
                                               double maxTimestepWidth,
                                               double endTime,
                                               LtsParameters lts)
    : vertexWeight(vertexWeight), cfl(cfl), maxTimestepWidth(maxTimestepWidth), endTime(endTime),
      lts(lts) {}

TimeSteppingParameters readTimeSteppingParameters(ParameterReader& baseReader) {
  auto reader = baseReader.readSubNode("discretization");
  const auto weightElement = reader.readWithDefault("vertexweightelement", 100);
  const auto weightDynamicRupture = reader.readWithDefault("vertexweightdynamicrupture", 100);
  const auto weightFreeSurfaceWithGravity =
      reader.readWithDefault("vertexweightfreesurfacewithgravity", 100);
  const double cfl = reader.readWithDefault("cfl", 0.5);
  const double maxTimestepWidth = reader.readWithDefault("fixtimestep", 5000.0);

  auto abortReader = baseReader.readSubNode("abortcriteria");
  const double endTime = reader.readWithDefault("endtime", 15.0);

  const LtsParameters ltsParameters = readLtsParameters(baseReader);

  reader.warnDeprecated({"ckmethod",
                         "dgfineout1d",
                         "fluxmethod",
                         "iterationcriterion",
                         "npoly",
                         "npolyrec",
                         "limitersecurityfactor",
                         "order",
                         "material",
                         "npolymap"});
  reader.warnUnknown();

  return TimeSteppingParameters({weightElement, weightDynamicRupture, weightFreeSurfaceWithGravity},
                                cfl,
                                maxTimestepWidth,
                                endTime,
                                ltsParameters);
}

} // namespace seissol::initializers::parameters
