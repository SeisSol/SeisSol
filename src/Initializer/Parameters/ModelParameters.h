// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_

#include "ParameterReader.h"

#include <string>

namespace seissol::initializer::parameters {

enum class ReflectionType { BothWaves = 1, BothWavesVelocity, Pwave, Swave };

struct ITMParameters {
  bool itmEnabled{false};
  double itmStartingTime{0.0};
  double itmDuration{0.0};
  double itmVelocityScalingFactor{1.0};
  ReflectionType itmReflectionType{ReflectionType::BothWaves};
};

enum class NumericalFlux { Godunov, Rusanov };

std::string fluxToString(NumericalFlux flux);

struct ModelParameters {
  bool hasBoundaryFile{false};
  bool plasticity{false};
  bool useCellHomogenizedMaterial{true};
  double freqCentral{};
  double freqRatio{1.0};
  double gravitationalAcceleration{};
  double tv{};
  std::string boundaryFileName;
  std::string materialFileName;
  std::vector<std::string> plasticityFileNames;
  ITMParameters itmParameters;
  NumericalFlux flux{NumericalFlux::Godunov};
  NumericalFlux fluxNearFault{NumericalFlux::Godunov};
};

ModelParameters readModelParameters(ParameterReader* baseReader);
ITMParameters readITMParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
