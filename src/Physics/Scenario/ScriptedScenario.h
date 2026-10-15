// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCENARIO_SCRIPTEDSCENARIO_H_
#define SEISSOL_SRC_PHYSICS_SCENARIO_SCRIPTEDSCENARIO_H_

#include "Common/ConfigRegistry.h"
#include "Physics/InitialField.h"
#include "Physics/ScriptField.h"

#include <initializer_list>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

// What the scenarios have in common that SeisSol writes as sderiv modules at setup: the constants
// they take from the material and the parameters are numbers in the source, which a field of the
// configuration then evaluates like a script of its own.

namespace seissol::physics::scenario {

/// `value` as an sderiv literal that reads back as the same double, the shortest one; in
/// parentheses if negative. A value that is not finite is a division that gives it.
std::string literal(double value);

/// `parts` one after the other.
std::string concat(std::initializer_list<std::string_view> parts);

/// The names of the components of the state of the configuration `config`.
std::vector<std::string> stateNamesOf(ConfigId config);

/// The names of the parameters of the material of the configuration `config`, as scripts read
/// them.
std::vector<std::string> parameterNamesOf(ConfigId config);

/// The fields of the sderiv module `source`, which SeisSol wrote for the scenario `name`, in the
/// configuration `config`: one per fused simulation if `perSimulation` (the module reads `sim`),
/// else one for all. `check`, if given, checks the material of every cell the fields are
/// evaluated in. A module that does not compile is an error that shows it.
std::vector<std::unique_ptr<InitialField>>
    scriptedFields(const std::string& source,
                   ConfigId config,
                   const std::string& name,
                   bool perSimulation,
                   const ScriptField::MaterialCheck& check = {});

} // namespace seissol::physics::scenario

#endif // SEISSOL_SRC_PHYSICS_SCENARIO_SCRIPTEDSCENARIO_H_
