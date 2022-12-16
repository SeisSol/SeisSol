#include "Initializer/InputAux.hpp"
#include "LtsParameters.h"
#include <utils/logger.h>

namespace seissol::initializers::time_stepping {

LtsParameters readLtsParametersFromYaml(std::shared_ptr<YAML::Node>& params) {
  using namespace seissol::initializers;

  const auto discretizationParams = (*params)["discretization"];
  const unsigned int rate = getWithDefault(discretizationParams, "clusteredlts", 1);
  const double wiggleFactorMinimum =
      getWithDefault(discretizationParams, "ltswigglefactormin", 1.0);
  const double wiggleFactorStepsize =
      getWithDefault(discretizationParams, "ltswigglefactorstepsize", 0.01);
  const bool wiggleFactorEnforceMaximumDifference =
      getWithDefault(discretizationParams, "ltswigglefactorenforcemaximumdifference", true);
  const unsigned int maxNumberOfClusters = getWithDefault(
      discretizationParams, "ltsmaxnumberofclusters", std::numeric_limits<int>::max() - 1);
  const double allowedRelativePerformanceLossAutoMerge =
      getWithDefault(discretizationParams, "ltsallowedrelativeperformancelossautomerge", 0.0);
  const double allowedPerformanceLossRatioAutoMerge = allowedRelativePerformanceLossAutoMerge + 1.0;
  return LtsParameters(rate,
                       wiggleFactorMinimum,
                       wiggleFactorStepsize,
                       wiggleFactorEnforceMaximumDifference,
                       maxNumberOfClusters,
                       allowedPerformanceLossRatioAutoMerge);
}

LtsParameters::LtsParameters(unsigned int rate,
                             double wiggleFactorMinimum,
                             double wiggleFactorStepsize,
                             bool wigleFactorEnforceMaximumDifference,
                             int maxNumberOfClusters,
                             double allowedPerformanceLossRatioAutoMerge)
    : rate(rate), wiggleFactorMinimum(wiggleFactorMinimum),
      wiggleFactorStepsize(wiggleFactorStepsize),
      wiggleFactorEnforceMaximumDifference(wigleFactorEnforceMaximumDifference),
      maxNumberOfClusters(maxNumberOfClusters),
      allowedPerformanceLossRatioAutoMerge(allowedPerformanceLossRatioAutoMerge) {
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

double LtsParameters::getWiggleFactorMinimum() const { return wiggleFactorMinimum; }

double LtsParameters::getWiggleFactorStepsize() const { return wiggleFactorStepsize; }

bool LtsParameters::getWiggleFactorEnforceMaximumDifference() const {
  return wiggleFactorEnforceMaximumDifference;
}

int LtsParameters::getMaxNumberOfClusters() const { return maxNumberOfClusters; }


double LtsParameters::getAllowedPerformanceLossRatioAutoMerge() const {
  return allowedPerformanceLossRatioAutoMerge;
}

} // namespace seissol::initializers::time_stepping
