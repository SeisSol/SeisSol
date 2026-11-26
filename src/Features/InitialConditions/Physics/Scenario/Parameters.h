// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_PARAMETERS_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_PARAMETERS_H_

#include "Core/Common/ConfigRegistry.h"
#include "WavePropagation/Physics/CommonDatastructures.h"

#include <Eigen/Dense>
#include <complex>
#include <functional>
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

namespace physics {

/// The plane wave operator of the material `material` of the configuration `config` for the wave
/// vector `k`, as the solver of the configuration lays the quantities out: a square matrix over
/// the quantities of the material, column by column. The plane waves take it from who sets them
/// up, since it depends on the solver (cf. solver::planeWaveOperator).
using PlaneWaveOperator = std::function<std::vector<std::complex<double>>(
    ConfigId config, const model::Material& material, const Eigen::Vector3d& k)>;

} // namespace physics

} // namespace seissol

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_PARAMETERS_H_
