// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_COMMON_CONFIGREGISTRY_H_
#define SEISSOL_SRC_CORE_COMMON_CONFIGREGISTRY_H_

#include "Core/Common/ConfigValue.h"
#include "Core/Common/Typedefs.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace seissol {

/**
 * @brief The index of a configuration among the ones built into the executable.
 *
 * It is also the variant of the configuration's kernels in runtime.h: the configuration of some
 * data, given by id, selects the kernels that compute on it.
 *
 * It depends on the build; outside of a run (e.g. in a file), a configuration is named by
 * `configName` instead.
 */
using ConfigId = std::uint32_t;

/**
 * @brief The quadrature rule of the dynamic rupture faces, the same for all configurations of the
 * executable.
 */
DRQuadRuleType builtDRQuadRule();

/// The number of configurations built into the executable.
std::size_t builtConfigCount();

/// The layout of a built configuration.
/// The solver that advances the cells of a built configuration, as the build chooses it for its
/// material.
SolverType configSolver(ConfigId id);

/// The value of a built configuration.
const ConfigValue& configValue(ConfigId id);

/// The id of a configuration, if it is built into the executable.
std::optional<ConfigId> findConfig(const ConfigValue& value);

/// The id of the configuration with the given canonical name, if it is built into the executable.
std::optional<ConfigId> findConfig(std::string_view name);

/// The names of the configurations built into the executable, in the order of their ids,
/// separated by commas.
std::string builtConfigNames();

/// The configuration that the environment variable SEISSOL_CONFIGURATION names, by its canonical
/// name, if it is set. Naming one that is not built into the executable is an error.
std::optional<ConfigId> environmentConfig();

/// The configuration of the cells for which the setup does not name one: the one that the
/// environment variable SEISSOL_CONFIGURATION names, or else the first one built.
ConfigId defaultConfig();

/// A description of all built configurations for humans: each name, followed by its fields and its
/// solver, and then what all of them share, the quadrature rule of the dynamic rupture faces.
std::string describeBuiltConfigs();

} // namespace seissol

#endif // SEISSOL_SRC_CORE_COMMON_CONFIGREGISTRY_H_
