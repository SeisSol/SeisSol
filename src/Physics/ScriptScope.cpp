// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Physics/ScriptScope.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Equations/Datastructures.h"
#include "Expr/Program.h"
#include "Expr/Rewrite.h"
#include "Model/ScriptParameters.h"
#include "Model/StateNames.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

namespace {

/// Whether a script keeps in `config` what depends on `name`: unless another configuration of the
/// run knows the name and `config` does not.
bool keeps(const std::vector<ConfigId>& configs,
           const KnowsName& knows,
           ConfigId config,
           const std::string& name) {
  if (knows(config, name)) {
    return true;
  }
  return std::none_of(configs.begin(), configs.end(), [&](ConfigId other) {
    return other != config && knows(other, name);
  });
}

} // namespace

std::vector<std::string> stateNamesOf(ConfigId config) {
  std::vector<std::string> names;
  dispatchConfig(config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    names = model::stateNames<model::MaterialOf<Cfg>>();
  });
  return names;
}

std::vector<std::string> quantityNamesOf(ConfigId config) {
  std::vector<std::string> names;
  dispatchConfig(config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    const auto& quantities = model::MaterialOf<Cfg>::Quantities;
    names.assign(quantities.begin(), quantities.end());
  });
  return names;
}

std::vector<std::string> parameterNamesOf(ConfigId config) {
  std::vector<std::string> names;
  dispatchConfig(config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    names = model::scriptParameterNames<model::MaterialOf<Cfg>>();
  });
  return names;
}

KnowsName knowsAmong(const std::vector<ConfigId>& configs,
                     const std::function<std::vector<std::string>(ConfigId)>& namesOf) {
  auto names = std::make_shared<std::map<ConfigId, std::set<std::string>>>();
  for (const auto config : configs) {
    const auto list = namesOf(config);
    (*names)[config] = std::set<std::string>(list.begin(), list.end());
  }
  return [names](ConfigId config, const std::string& name) {
    const auto found = names->find(config);
    return found != names->end() && found->second.count(name) > 0;
  };
}

ScriptScope::ScriptScope(std::vector<ConfigId> configs, KnowsName provides, KnowsName takes)
    : configs_(std::move(configs)), provides_(std::move(provides)), takes_(std::move(takes)) {}

bool ScriptScope::keepsInput(ConfigId config, const std::string& name) const {
  return keeps(configs_, provides_, config, name);
}

bool ScriptScope::keepsOutput(ConfigId config, const std::string& name) const {
  return keeps(configs_, takes_, config, name);
}

expr::Program ScriptScope::restrict(const expr::Program& program,
                                    ConfigId config,
                                    const std::string& name,
                                    bool log) const {
  expr::Restriction restriction;
  auto restricted = expr::restrictProgram(
      program,
      [&](const std::string& input) { return keepsInput(config, input); },
      [&](const std::string& output) { return keepsOutput(config, output); },
      &restriction);
  if (log && (!restriction.unwanted.empty() || !restriction.unavailable.empty())) {
    std::ostringstream message;
    for (const auto& output : restriction.unwanted) {
      message << " " << output << " (of another material)";
    }
    for (const auto& [output, input] : restriction.unavailable) {
      message << " " << output << " (reads " << input << ", of another material)";
    }
    logInfo() << "The script" << name << "leaves out in the configuration"
              << configName(configValue(config)) << ":" << message.str();
  }
  return restricted;
}

ProgramFilter ScriptScope::filter(ConfigId config, std::string name, bool log) const {
  return [scope = *this, config, name = std::move(name), log](const expr::Program& program) {
    return scope.restrict(program, config, name, log);
  };
}

std::unique_ptr<reader::scripting::DataReader> ScriptScope::reader(
    const std::string& path, const std::vector<std::string>& defaultInArgs, ConfigId config) const {
  auto reader = reader::scripting::buildReader(path, defaultInArgs);
  const auto* compiled = dynamic_cast<const reader::scripting::CompiledReader*>(reader.get());
  if (compiled == nullptr) {
    return reader;
  }
  auto restricted = restrict(compiled->program(), config, path);
  if (restricted.outputs().size() == compiled->program().outputs().size()) {
    return reader;
  }
  // a reference of the whole script would read what the configuration does not have
  return std::make_unique<reader::scripting::CompiledReader>(std::move(restricted), nullptr);
}

void ScriptScope::checkKept(const expr::Program& program, const std::string& name) const {
  std::set<std::string> kept;
  for (const auto config : configs_) {
    const auto restricted = expr::restrictProgram(
        program,
        [&](const std::string& input) { return keepsInput(config, input); },
        [&](const std::string& output) { return keepsOutput(config, output); });
    for (const auto& output : restricted.outputs()) {
      kept.insert(output.name);
    }
  }
  for (const auto& output : program.outputs()) {
    if (kept.count(output.name) == 0) {
      logError() << "The script" << name << "gives" << output.name
                 << "in none of the configurations of the run: it reads or gives names of "
                    "several materials, which no configuration has all of.";
    }
  }
}

ScriptScope fieldScope(const std::vector<ConfigId>& configs) {
  return {configs, knowsAmong(configs, parameterNamesOf), knowsAmong(configs, stateNamesOf)};
}

ScriptScope ghostStateScope(const std::vector<ConfigId>& configs) {
  return {configs,
          knowsAmong(configs,
                     [](ConfigId config) {
                       auto names = quantityNamesOf(config);
                       const auto parameters = parameterNamesOf(config);
                       names.insert(names.end(), parameters.begin(), parameters.end());
                       return names;
                     }),
          knowsAmong(configs, quantityNamesOf)};
}

} // namespace seissol::physics
