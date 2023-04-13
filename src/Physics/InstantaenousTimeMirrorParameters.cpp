#include "Initializer/InputAux.hpp"
#include "InstantaneousTimeMirrorParameters.h"
#include <utils/logger.h>

namespace seissol::initializers::ITM {
ITMParameters readITMParametersFromYaml(std::shared_ptr<YAML::Node>& params) {
  using namespace seissol::initializers;

  const auto equationsParams = (*params)["equations"];
  const double ITMTime = getWithDefault(equationsParams, "itmtime", 1.0);

  return ITMParameters(ITMTime);
}

ITMParameters::ITMParameters(double ITMTime) : ITMTime(ITMTime) {
  if (ITMTime < 0.0) {
    logError() << "ITM Time is less than zero. It should be positive!\n";
  }
}

double ITMParameters::getITMTime() const{return ITMTime;}

} // namespace seissol::initializers::ITM