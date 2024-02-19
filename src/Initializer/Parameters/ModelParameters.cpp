#include "ModelParameters.h"

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
DamagedElasticParameters readDamagedElasticParameters(ParameterReader* baseReader){
    auto* reader = baseReader -> readSubNode("equations");
    const auto epsInitxx = reader->readWithDefault<real>("epsinitxx", 0);
    const auto epsInityy = reader->readWithDefault<real>("epsinityy", 0);
    const auto epsInitzz = reader->readWithDefault<real>("epsinitzz", 0);
    const auto epsInitxy = reader->readWithDefault<real>("epsinitxy", 10);
    const auto epsInityz = reader->readWithDefault<real>("epsinityz", 0);
    const auto epsInitzx = reader->readWithDefault<real>("epsinitzx", 0);
    const auto beta_alpha = reader->readWithDefault<real>("betaalpha", 0);
    const auto aB0 = reader->readWithDefault<real>("aB0", 0);
    const auto aB1 = reader->readWithDefault<real>("aB1", 0);
    const auto aB2 = reader->readWithDefault<real>("aB2", 0);
    const auto aB3 = reader->readWithDefault<real>("aB3", 0);
    const auto scalingvalue = reader->readWithDefault<real>("scalingvalue", 0);

    if(!isModelDamagedElastic){
      reader->markUnused(
          {"epsinitxx", "epsinityy", "epsinitzz", "epsinitxy", "epsinityz", "epsinitzx", "betaalpha", "aB0", "aB1", "aB2", "aB3", "scalingvalue"});
    }
    return DamagedElasticParameters{
      epsInitxx, epsInityy, epsInitzz, epsInitxy, epsInityz, epsInitzx, beta_alpha, aB0, aB1, aB2, aB3, scalingvalue};
}

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

  const ITMParameters itmParameters = readITMParameters(baseReader);
  const DamagedElasticParameters damagedElasticParameters = readDamagedElasticParameters(baseReader);

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  return ModelParameters{hasBoundaryFile,
                         plasticity,
                         useCellHomogenizedMaterial,
                         freqCentral,
                         freqRatio,
                         gravitationalAcceleration,
                         tv,
                         boundaryFileName,
                         materialFileName,
                         itmParameters,
                         damagedElasticParameters
  };
}
} // namespace seissol::initializer::parameters
