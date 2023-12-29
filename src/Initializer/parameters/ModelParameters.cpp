#include "ModelParameters.h"

namespace seissol::initializers::parameters {

ModelParameters readModelParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("equations");

  const std::string boundaryFileName = reader->readWithDefault("boundaryfileName", std::string(""));
  const std::string materialFileName =
      reader->readOrFail<std::string>("materialfilename", "No material file given.");
  const bool hasBoundaryFile = boundaryFileName != "";

  const bool plasticity = reader->readWithDefault("plasticity", false);
  const bool useCellHomogenizedMaterial =
      reader->readWithDefault("usecellhomogenizedmaterial", true);

  const double gravitationalAcceleration =
      reader->readWithDefault("gravitationalacceleration", 9.81);
  const double tv = reader->readWithDefault("tv", 0.1);

  const double freqCentral = reader->readIfRequired<double>("freqcentral", isModelViscoelastic());
  const double freqRatio = reader->readIfRequired<double>("freqratio", isModelViscoelastic());
  if constexpr (isModelViscoelastic()) {
    if (freqRatio <= 0) {
      logError()
          << "The freqratio parameter must be positive---but that is currently not the case.";
    }
  }

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  return ModelParameters{hasBoundaryFile,
                         plasticity,
                         useCellHomogenizedMaterial,
                         freqCentral,
                         freqRatio,
                         gravitationalAcceleration,
                         tv,
                         boundaryFileName,
                         materialFileName};
}
} // namespace seissol::initializers::parameters
