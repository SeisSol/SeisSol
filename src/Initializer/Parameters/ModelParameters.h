// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_

#include <string>

#include "ParameterReader.h"

namespace seissol::initializer::parameters {

enum class ReflectionType { BothWaves = 1, BothWavesVelocity, Pwave, Swave };

struct ITMParameters {
  bool itmEnabled;
  double itmStartingTime;
  double itmDuration;
  double itmVelocityScalingFactor;
  ReflectionType itmReflectionType;
};

enum class NumericalFlux { Godunov, Rusanov };

std::string fluxToString(NumericalFlux flux);

struct ModelParameters {
  bool hasBoundaryFile;
  bool plasticity;
  bool useCellHomogenizedMaterial;
  double freqCentral;
  double freqRatio;
  double gravitationalAcceleration;
  double tv;
  std::string boundaryFileName;
  std::string materialFileName;
  ITMParameters itmParameters;
  NumericalFlux flux;
  NumericalFlux fluxNearFault;
};

ModelParameters readModelParameters(ParameterReader* baseReader);
ITMParameters readITMParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
