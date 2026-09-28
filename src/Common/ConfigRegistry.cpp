// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigRegistry.h"

#include "Common/ConfigLayout.h"
#include "Common/ConfigValue.h"

#include <cstddef>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>

namespace seissol {

std::size_t builtConfigCount() { return builtConfigLayouts().size(); }

const ConfigLayout& configLayout(ConfigId id) { return builtConfigLayouts().at(id); }

const ConfigValue& configValue(ConfigId id) { return configLayout(id).config; }

std::optional<ConfigId> findConfig(const ConfigValue& value) {
  const auto& layouts = builtConfigLayouts();
  for (std::size_t id = 0; id < layouts.size(); ++id) {
    if (layouts[id].config == value) {
      return static_cast<ConfigId>(id);
    }
  }
  return std::nullopt;
}

std::optional<ConfigId> findConfig(std::string_view name) {
  const auto value = parseConfigName(name);
  if (!value.has_value()) {
    return std::nullopt;
  }
  return findConfig(value.value());
}

ConfigId defaultConfig() { return 0; }

std::string describeBuiltConfigs() {
  std::ostringstream stream;
  for (std::size_t id = 0; id < builtConfigCount(); ++id) {
    const auto& value = configValue(static_cast<ConfigId>(id));
    stream << configName(value) << "\n";
    std::istringstream fields(describeConfig(value));
    for (std::string field; std::getline(fields, field);) {
      stream << "  " << field << "\n";
    }
  }
  return stream.str();
}

} // namespace seissol
