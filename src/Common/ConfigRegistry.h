// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONFIGREGISTRY_H_
#define SEISSOL_SRC_COMMON_CONFIGREGISTRY_H_

#include "Common/ConfigLayout.h"
#include "Common/ConfigValue.h"

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
 * @brief The layouts of all configurations built into the executable, indexed by their id.
 *
 * Provided by the code that is generated for the configurations; everything else goes through
 * the functions below.
 */
const std::vector<ConfigLayout>& builtConfigLayouts();

/// The number of configurations built into the executable.
std::size_t builtConfigCount();

/// The layout of a built configuration.
const ConfigLayout& configLayout(ConfigId id);

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

/// A description of all built configurations for humans: each name, followed by its fields.
std::string describeBuiltConfigs();

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONFIGREGISTRY_H_
