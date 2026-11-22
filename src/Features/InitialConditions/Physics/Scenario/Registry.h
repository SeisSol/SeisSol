// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_REGISTRY_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_REGISTRY_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Memory/Typedefs.h"
#include "Features/InitialConditions/InitializationParameters.h"
#include "Features/InitialConditions/Physics/InitialField.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

namespace seissol::physics::scenario {

/**
 * What a scenario with a time mirror takes of the settings of the mirror: whether the run has
 * one, when it starts and how long it lasts, and the factor it scales the velocities by.
 */
struct MirrorSettings {
  bool enabled{false};
  double startingTime{0};
  double duration{0};
  double velocityScalingFactor{1};
};

/**
 * Everything a scenario needs to be constructed. Bundled so that adding a scenario with new
 * requirements does not change the signature of every entry in the registry.
 */
struct Input {
  const initializer::parameters::InitializationParameters& initialization;
  // the settings of the time mirror, which the scenario with ITM passes through
  MirrorSettings timeMirror;
  // the number of simulations of the run
  std::size_t simulationCount{1};
  const CellMaterialData& materialData;
  GravitationSetup gravitation;
  // the configuration the scenario is set up for: its material and its fused simulations
  ConfigId config{};
  // the simulation of the run to set the scenario up for alone, e.g. with a material of its own;
  // all of them if not given
  std::optional<std::size_t> simulation;
};

/**
 * Whether a scenario can be constructed in this configuration, and if not, why.
 */
struct Availability {
  bool available{true};
  std::string_view reason;
};

/**
 * The name the scenario is reported under.
 */
std::string_view name(initializer::parameters::InitializationType type);

/**
 * Whether the scenario is known and defined for the material of the configuration `config`.
 */
Availability availability(initializer::parameters::InitializationType type, ConfigId config);

/**
 * Whether the configured scenario can serve an analytical boundary condition, which needs the
 * solution at arbitrary times rather than only at t = 0.
 */
Availability analyticalBoundaryAvailability(initializer::parameters::InitializationType type,
                                            ConfigId config);

/**
 * The acoustic travelling wave with ITM as the parameter file sets it up. The wave passes the
 * mirror only if ITM is enabled; otherwise, it travels on unchanged, whatever the mirror settings.
 */
AcousticTravellingWaveParametersITM acousticTravellingWaveParameters(
    const initializer::parameters::InitializationParameters& initialization,
    const MirrorSettings& mirror);

/**
 * Constructs the scenario: one instance per simulation of the run, or a single one for all of them
 * where they do not differ; for `Input::simulation`, a single one in any case. Only call after
 * availability() reported it as available.
 */
std::vector<std::unique_ptr<InitialField>> build(initializer::parameters::InitializationType type,
                                                 const Input& input);

} // namespace seissol::physics::scenario

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_PHYSICS_SCENARIO_REGISTRY_H_
