// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Physics/Scenario/ScriptedScenario.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Physics/InitialField.h"
#include "Physics/ScriptField.h"

#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics::scenario {

std::string literal(double value) {
  // sderiv has no literal for a value that is not finite; a division gives it
  if (std::isnan(value)) {
    return "(0.0 / 0.0)";
  }
  if (std::isinf(value)) {
    return value > 0 ? "(1.0 / 0.0)" : "(-1.0 / 0.0)";
  }
  std::array<char, 32> buffer{};
  const auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), std::abs(value));
  std::string digits(buffer.data(), result.ptr);
  if (!std::signbit(value)) {
    return digits;
  }
  std::string negative = "(-";
  negative.append(digits).append(")");
  return negative;
}

std::string concat(std::initializer_list<std::string_view> parts) {
  std::string result;
  for (const auto& part : parts) {
    result.append(part);
  }
  return result;
}

std::vector<std::unique_ptr<InitialField>> scriptedFields(const std::string& source,
                                                          ConfigId config,
                                                          const std::string& name,
                                                          bool perSimulation,
                                                          const ScriptField::MaterialCheck& check) {
  std::optional<expr::Program> program;
  try {
    program = expr::compileSderivModule(source);
  } catch (const expr::SderivError& error) {
    logError() << "The script of" << name << "does not compile:" << error.what() << "\n" << source;
  }

  const auto names = stateNamesOf(config);
  const auto parameters = parameterNamesOf(config);
  std::vector<std::unique_ptr<InitialField>> fields;
  const std::size_t count = perSimulation ? configValue(config).numSimulations : 1;
  for (std::size_t simulation = 0; simulation < count; ++simulation) {
    auto field = std::make_unique<ScriptField>(*program, name, names, parameters, simulation);
    if (check) {
      field->checkMaterial(check);
    }
    fields.push_back(std::move(field));
  }
  return fields;
}

} // namespace seissol::physics::scenario
