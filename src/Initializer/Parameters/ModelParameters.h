// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_

#include "ParameterReader.h"

#include <array>
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

/// Parameters of a continuum damage-breakage rheology that belong to the
/// model rather than to a point in the mesh, so they come from here rather
/// than from easi. The defaults are the configuration the published
/// implementation runs: no breakage, no healing, and a granular branch that
/// is never reached.
///
/// What those defaults mean for a run: the damage grows to the critical damage
/// of the cell's strain direction and stops there, because that is where the
/// solid branch stops describing anything and there is no breakage to take
/// over. The state past that point is not modeled rather than modeled badly,
/// which is the side to err on -- the medium has no waves there.
struct DamageParameters {
  double breakageRate{0.0};
  double healingRate{0.0};
  /// Width of the smoothed step from damage to breakage. It divides, so it
  /// must not be zero even where breakage is switched off.
  double betaAlpha{1.0};
  /// Coefficients of the granular branch, aB0 to aB3.
  std::array<double, 4> granular{};
};

enum class NumericalFlux { Godunov, Rusanov };

std::string fluxToString(NumericalFlux flux);

struct ModelParameters {
  bool hasBoundaryFile{false};
  bool plasticity{false};
  bool plasticityPointwise{true};
  std::unordered_set<int> plasticityDisabledGroups;
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
  DamageParameters damageParameters;
};

ModelParameters readModelParameters(ParameterReader* baseReader);
ITMParameters readITMParameters(ParameterReader* baseReader);
DamageParameters readDamageParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_MODELPARAMETERS_H_
