#include "Initializer/InputAux.hpp"
#include "LtsParameters.h"
#include <utils/logger.h>

seissol::time_stepping::LtsParameters
    seissol::time_stepping::readLtsParametersFromYaml(std::shared_ptr<YAML::Node>& params) {
  using namespace seissol::initializers;

  const auto discretizationParams = (*params)["discretization"];
  unsigned int const rate = getWithDefault(discretizationParams, "clusteredlts", 1);
  const double wiggleFactorMinimum =
      getWithDefault(discretizationParams, "ltswigglefactormin", 1.0);
  const double wiggleFactorStepsize =
      getWithDefault(discretizationParams, "ltswigglefactorstepsize", 0.01);
  return seissol::time_stepping::LtsParameters(rate, wiggleFactorMinimum, wiggleFactorStepsize);
}

seissol::time_stepping::LtsParameters::LtsParameters(unsigned int rate,
                                                     double wiggleFactorMinimum,
                                                     double wiggleFactorStepsize)
    : rate(rate), wiggleFactorMinimum(wiggleFactorMinimum),
      wiggleFactorStepsize(wiggleFactorStepsize) {
  const bool isWiggleFactorValid = wiggleFactorMinimum <= 1.0 && wiggleFactorMinimum > (1.0 / rate);
  if (!isWiggleFactorValid) {
    logError() << "Minimal wiggle factor of " << wiggleFactorMinimum << "is not valid for rate"
               << rate;
  }
}

bool seissol::time_stepping::LtsParameters::isWiggleFactorUsed() const {
  return wiggleFactorMinimum < 1.0;
}
unsigned int seissol::time_stepping::LtsParameters::getRate() const { return rate; }
double seissol::time_stepping::LtsParameters::getWiggleFactorMinimum() const {
  return wiggleFactorMinimum;
}
double seissol::time_stepping::LtsParameters::getWiggleFactorStepsize() const {
  return wiggleFactorStepsize;
}
