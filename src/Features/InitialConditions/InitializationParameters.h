// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_INITIALIZATIONPARAMETERS_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_INITIALIZATIONPARAMETERS_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

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
  InitializationType type{InitializationType::Zero};
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  // per quantity of the material of the configuration of the run
  std::vector<double> ampField;
  double magnitude{};
  double width{};
  double k{};
  // the easi file of every simulation of the run, the one of the first if not given
  std::vector<std::string> filenames;
  bool hasTime{false};
  bool avoidIC{false};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_INITIALIZATIONPARAMETERS_H_
