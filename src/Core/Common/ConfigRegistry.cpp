// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigRegistry.h"

#include "Config.h"
#include "Core/Common/ConfigValue.h"
#include "Core/Common/Typedefs.h"

#include <cstddef>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol {

namespace {

/// A configuration built into the executable, as Config.h lists them: in the order of their ids.
struct BuiltConfig {
  ConfigValue config;
  SolverType solver{};
};

const std::vector<BuiltConfig>& builtConfigs() {
#define SEISSOL_BUILT_CONFIG(Cfg) BuiltConfig{Cfg::Value, Cfg::Solver},
  static const std::vector<BuiltConfig> Configs{SEISSOL_FOR_EACH_CONFIG(SEISSOL_BUILT_CONFIG)};
#undef SEISSOL_BUILT_CONFIG
  return Configs;
}

} // namespace

DRQuadRuleType builtDRQuadRule() { return DRQuadRule; }

std::size_t builtConfigCount() { return builtConfigs().size(); }

const ConfigValue& configValue(ConfigId id) { return builtConfigs().at(id).config; }

SolverType configSolver(ConfigId id) { return builtConfigs().at(id).solver; }

std::optional<ConfigId> findConfig(const ConfigValue& value) {
  const auto& configs = builtConfigs();
  for (std::size_t id = 0; id < configs.size(); ++id) {
    if (configs[id].config == value) {
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
    const auto& built = builtConfigs()[id];
    stream << configName(built.config) << "\n";
    std::istringstream fields(describeConfig(built.config));
    for (std::string field; std::getline(fields, field);) {
      stream << "  " << field << "\n";
    }
    stream << "  Solver: " << solverTypeName(built.solver) << "\n";
  }
  stream << "Dynamic rupture quadrature rule (all configurations): "
         << drQuadRuleName(builtDRQuadRule()) << "\n";
  return stream.str();
}

} // namespace seissol
