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
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>

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

std::string builtConfigNames() {
  std::string names;
  for (std::size_t id = 0; id < builtConfigCount(); ++id) {
    names += (id > 0 ? ", " : "") + configName(configValue(static_cast<ConfigId>(id)));
  }
  return names;
}

std::optional<ConfigId> environmentConfig() {
  auto name = utils::Env("SEISSOL_").getOptional<std::string>("CONFIGURATION");
  if (!name.has_value()) {
    return std::nullopt;
  }
  // spelled as in a parameter file
  utils::StringUtils::trim(name.value());
  utils::StringUtils::toLower(name.value());
  const auto config = findConfig(name.value());
  if (!config.has_value()) {
    logError() << "The configuration" << name.value()
               << "of SEISSOL_CONFIGURATION is not built into this executable. It has:"
               << builtConfigNames();
  }
  return config;
}

ConfigId defaultConfig() { return environmentConfig().value_or(0); }

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
