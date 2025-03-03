// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ModelParameters.h"
#include <Equations/Datastructures.h>
#include <Initializer/Parameters/ParameterReader.h>
#include <Model/CommonDatastructures.h>
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

  constexpr auto IsViscoelastic = model::MaterialT::Type == model::MaterialType::Viscoelastic;

  const auto freqCentral = reader->readIfRequired<double>("freqcentral", IsViscoelastic);
  const auto freqRatio = reader->readIfRequired<double>("freqratio", IsViscoelastic);
  if constexpr (IsViscoelastic) {
    if (freqRatio <= 0) {
      logError() << "The freqratio parameter must be positive; but that is currently not the case.";
    }
  }

  const ITMParameters itmParameters = readITMParameters(baseReader);

  reader->warnDeprecated({"adjoint", "adjfilename", "anisotropy"});

  const auto flux =
      reader->readWithDefaultStringEnum<NumericalFlux>("numflux",
                                                       "godunov",
                                                       {
                                                           {"godunov", NumericalFlux::Godunov},
                                                           {"rusanov", NumericalFlux::Rusanov},
                                                       });

  const auto fluxNearFault =
      reader->readWithDefaultStringEnum<NumericalFlux>("numfluxnearfault",
                                                       "godunov",
                                                       {
                                                           {"godunov", NumericalFlux::Godunov},
                                                           {"rusanov", NumericalFlux::Rusanov},
                                                       });

  return ModelParameters{hasBoundaryFile,
                         plasticity,
                         useCellHomogenizedMaterial,
                         freqCentral,
                         freqRatio,
                         gravitationalAcceleration,
                         tv,
                         boundaryFileName.value_or(""),
                         materialFileName,
                         itmParameters,
                         flux,
                         fluxNearFault};
}

std::string fluxToString(NumericalFlux flux) {
  if (flux == NumericalFlux::Godunov) {
    return "Godunov flux";
  }
  if (flux == NumericalFlux::Rusanov) {
    return "Rusanov flux";
  }
  return "(unknown flux)";
}

} // namespace seissol::initializer::parameters
