#include "Initializer/InputAux.hpp"
#include "LtsParameters.h"
#include <utils/logger.h>

namespace seissol::initializers::time_stepping {

LtsParameters readLtsParametersFromYaml(std::shared_ptr<YAML::Node>& params) {
  using namespace seissol::initializers;

  const auto discretizationParams = (*params)["discretization"];
  unsigned int const rate = getWithDefault(discretizationParams, "clusteredlts", 1);
  const double wiggleFactorMinimum =
      getWithDefault(discretizationParams, "ltswigglefactormin", 1.0);
  const double wiggleFactorStepsize =
      getWithDefault(discretizationParams, "ltswigglefactorstepsize", 0.01);
  const bool wiggleFactorEnforceMaximumDifference =
      getWithDefault(discretizationParams, "ltswigglefactorenforcemaximumdifference", true);
  unsigned int const maxNumberOfClusters = getWithDefault(
      discretizationParams, "ltsmaxnumberofclusters", std::numeric_limits<int>::max() - 1);
  return LtsParameters(rate,
                       wiggleFactorMinimum,
                       wiggleFactorStepsize,
                       wiggleFactorEnforceMaximumDifference,
                       maxNumberOfClusters);
}

LtsParameters::LtsParameters(unsigned int rate,
                             double wiggleFactorMinimum,
                             double wiggleFactorStepsize,
                             bool wigleFactorEnforceMaximumDifference,
                             int maxNumberOfClusters)
    : rate(rate), wiggleFactorMinimum(wiggleFactorMinimum),
      wiggleFactorStepsize(wiggleFactorStepsize),
      wiggleFactorEnforceMaximumDifference(wigleFactorEnforceMaximumDifference),
      maxNumberOfClusters(maxNumberOfClusters) {
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
}

bool LtsParameters::isWiggleFactorUsed() const { return wiggleFactorMinimum < 1.0; }

unsigned int LtsParameters::getRate() const { return rate; }

double LtsParameters::getWiggleFactorMinimum() const { return wiggleFactorMinimum; }

double LtsParameters::getWiggleFactorStepsize() const { return wiggleFactorStepsize; }

bool LtsParameters::getWiggleFactorEnforceMaximumDifference() const {
  return wiggleFactorEnforceMaximumDifference;
}

int LtsParameters::getMaxNumberOfClusters() const { return maxNumberOfClusters; }

} // namespace seissol::initializers::time_stepping
