#include "ModelParameters.h"
#include <Initializer/Parameters/ParameterReader.h>
#include <utils/logger.h>

namespace seissol::initializer::parameters {

ITMParameters readITMParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("equations");
  const auto itmEnabled = reader->readWithDefault<bool>("itmenable", false);
  const auto itmStartingTime = reader->readWithDefault<double>("itmstartingtime", 0.0);
  const auto itmDuration = reader->readWithDefault<double>("itmtime", 0.0);
  const auto itmVelocityScalingFactor =
      reader->readWithDefault<double>("itmvelocityscalingfactor", 1.0);
  const auto reflectionType =
      reader->readWithDefaultEnum<ReflectionType>("itmreflectiontype",
                                                  ReflectionType::BothWaves,
                                                  {ReflectionType::BothWaves,
                                                   ReflectionType::BothWavesVelocity,
                                                   ReflectionType::Pwave,
                                                   ReflectionType::Swave});
  if (itmEnabled) {
    if (itmDuration <= 0.0) {
      logError() << "ITM Time is not positive. It should be positive!";
    }
    if (itmVelocityScalingFactor < 0.0) {
      logError() << "ITM Velocity Scaling Factor is less than zero. It should be positive!";
    }
    if (itmStartingTime < 0.0) {
      logError() << "ITM Starting Time can not be less than zero";
    }
  } else {
    reader->markUnused(
        {"itmstartingtime", "itmtime", "itmvelocityscalingfactor", "itmreflectiontype"});
  }
  return ITMParameters{
      itmEnabled, itmStartingTime, itmDuration, itmVelocityScalingFactor, reflectionType};
}

ModelParameters readModelParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("equations");

  const auto boundaryFileName = reader->readPath("boundaryfileName");
  const std::string materialFileName =
      reader->readPathOrFail("materialfilename", "No material file given.");
  const bool hasBoundaryFile = !boundaryFileName.value_or("").empty();

  const bool plasticity = reader->readWithDefault("plasticity", false);
  const bool useCellHomogenizedMaterial =
      reader->readWithDefault("usecellhomogenizedmaterial", true);

  const double gravitationalAcceleration =
      reader->readWithDefault("gravitationalacceleration", 9.81);
  const double tv = reader->readWithDefault("tv", 0.1);

  const auto freqCentral = reader->readIfRequired<double>("freqcentral", isModelViscoelastic());
  const auto freqRatio = reader->readIfRequired<double>("freqratio", isModelViscoelastic());
  if constexpr (isModelViscoelastic()) {
    if (freqRatio <= 0) {
      logError()
          << "The freqratio parameter must be positive---but that is currently not the case.";
    }
  }

  const ITMParameters itmParameters = readITMParameters(baseReader);

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  return ModelParameters{hasBoundaryFile,
                         plasticity,
                         useCellHomogenizedMaterial,
                         freqCentral,
                         freqRatio,
                         gravitationalAcceleration,
                         tv,
                         boundaryFileName.value_or(""),
                         materialFileName,
                         itmParameters};
}
} // namespace seissol::initializer::parameters
