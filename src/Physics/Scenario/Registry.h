// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCENARIO_REGISTRY_H_
#define SEISSOL_SRC_PHYSICS_SCENARIO_REGISTRY_H_

#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Physics/InitialField.h"

#include <memory>
#include <string_view>
#include <vector>

namespace seissol::physics::scenario {

/**
 * Everything a scenario needs to be constructed. Bundled so that adding a scenario with new
 * requirements does not change the signature of every entry in the registry.
 */
struct Input {
  const initializer::parameters::SeisSolParameters& parameters;
  const CellMaterialData& materialData;
  GravitationSetup gravitation;
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
 * Whether the scenario is known and defined for the material this binary was built for.
 */
Availability availability(initializer::parameters::InitializationType type);

/**
 * Whether the configured scenario can serve an analytical boundary condition, which needs the
 * solution at arbitrary times rather than only at t = 0.
 */
Availability analyticalBoundaryAvailability(initializer::parameters::InitializationType type);

/**
 * Constructs the scenario, one instance per fused simulation. Only call after availability()
 * reported it as available.
 */
std::vector<std::unique_ptr<InitialField>> build(initializer::parameters::InitializationType type,
                                                 const Input& input);

} // namespace seissol::physics::scenario

#endif // SEISSOL_SRC_PHYSICS_SCENARIO_REGISTRY_H_
