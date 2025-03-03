// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_INITIALIZATIONPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_INITIALIZATIONPARAMETERS_H_

#include <Eigen/Dense>

#include "Equations/Datastructures.h"
#include "Initializer/InputAux.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

enum class InitializationType : int {
  Zero,
  Planarwave,
  SuperimposedPlanarwave,
  Travelling,
  AcousticTravellingWithITM,
  Scholte,
  Snell,
  Ocean0,
  Ocean1,
  Ocean2,
  PressureInjection,
  Easi
};

struct InitializationParameters {
  InitializationType type;
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  Eigen::Vector<double, seissol::model::MaterialT::NumQuantities> ampField;
  double magnitude;
  double width;
  double k;
  std::string filename;
  bool hasTime;
};

InitializationParameters readInitializationParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_INITIALIZATIONPARAMETERS_H_
