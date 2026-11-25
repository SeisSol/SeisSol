// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_PARAMETERS_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_PARAMETERS_H_

#include <Eigen/Dense>
#include <complex>
#include <vector>

namespace seissol {

/// a plane wave: its origin, its wave vector, and the quantities it carries with their amplitudes
struct TravellingWaveParameters {
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  std::vector<int> varField;
  std::vector<std::complex<double>> ampField;
};

/// an acoustic travelling wave through an instantaneous time mirror
struct AcousticTravellingWaveParametersITM {
  double k{};
  double itmStartingTime{};
  double itmDuration{};
  double itmVelocityScalingFactor{};
};

} // namespace seissol

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_PARAMETERS_H_
