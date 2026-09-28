// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigValue.h"

#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Model/MaterialType.h"

#include <array>
#include <cstddef>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace seissol {

namespace {

template <typename T>
using NameEntry = std::pair<T, std::string_view>;

constexpr std::array<NameEntry<model::MaterialType>, 7> MaterialNames{{
    {model::MaterialType::Solid, "solid"},
    {model::MaterialType::Acoustic, "acoustic"},
    {model::MaterialType::Elastic, "elastic"},
    {model::MaterialType::Viscoelastic, "viscoelastic"},
    {model::MaterialType::Viscoacoustic, "viscoacoustic"},
    {model::MaterialType::Anisotropic, "anisotropic"},
    {model::MaterialType::Poroelastic, "poroelastic"},
}};

constexpr std::array<NameEntry<SolverType>, 3> SolverNames{{
    {SolverType::LinearCK, "linearck"},
    {SolverType::LinearCKAnelastic, "linearckanelastic"},
    {SolverType::STP, "stp"},
}};

constexpr std::array<NameEntry<DRQuadRuleType>, 3> DRQuadRuleNames{{
    {DRQuadRuleType::Stroud, "stroud"},
    {DRQuadRuleType::Dunavant, "dunavant"},
    {DRQuadRuleType::WitherdenVincent, "witherdenvincent"},
}};

constexpr std::array<RealType, 2> Precisions{RealType::F32, RealType::F64};

template <typename T, std::size_t N>
std::string_view nameOf(const std::array<NameEntry<T>, N>& table, T value) {
  for (const auto& [entry, name] : table) {
    if (entry == value) {
      return name;
    }
  }
  return "unknown";
}

template <typename T, std::size_t N>
std::optional<T> entryOf(const std::array<NameEntry<T>, N>& table, std::string_view name) {
  for (const auto& [entry, entryName] : table) {
    if (entryName == name) {
      return entry;
    }
  }
  return std::nullopt;
}

std::optional<RealType> precisionOf(std::string_view name) {
  for (const auto precision : Precisions) {
    if (stringRealType(precision) == name) {
      return precision;
    }
  }
  return std::nullopt;
}

/// A number after a one-letter prefix, e.g. `o6`; nothing but digits may follow the prefix.
std::optional<std::size_t> prefixedNumber(std::string_view token, char prefix) {
  if (token.size() < 2 || token.front() != prefix) {
    return std::nullopt;
  }
  std::size_t number = 0;
  for (const char digit : token.substr(1)) {
    if (digit < '0' || digit > '9') {
      return std::nullopt;
    }
    // a number too large to be represented does not survive the comparison with the canonical
    // name in `parseConfigName`
    number = 10 * number + static_cast<std::size_t>(digit - '0');
  }
  return number;
}

std::vector<std::string_view> split(std::string_view text, char separator) {
  std::vector<std::string_view> parts;
  std::size_t start = 0;
  while (true) {
    const auto position = text.find(separator, start);
    if (position == std::string_view::npos) {
      parts.push_back(text.substr(start));
      return parts;
    }
    parts.push_back(text.substr(start, position - start));
    start = position + 1;
  }
}

} // namespace

std::string_view materialTypeName(model::MaterialType type) { return nameOf(MaterialNames, type); }

std::string_view solverTypeName(SolverType type) { return nameOf(SolverNames, type); }

std::string_view drQuadRuleName(DRQuadRuleType rule) { return nameOf(DRQuadRuleNames, rule); }

std::string configName(const ConfigValue& value) {
  std::string name;
  name += materialTypeName(value.materialType);
  name += '-';
  name += solverTypeName(value.solver);
  if (value.relaxationMechanisms > 0) {
    name += "-m" + std::to_string(value.relaxationMechanisms);
  }
  name += "-o" + std::to_string(value.convergenceOrder);
  name += '-';
  name += stringRealType(value.precision);
  name += '-';
  name += drQuadRuleName(value.drQuadRule);
  if (value.numSimulations != 1) {
    name += "-s" + std::to_string(value.numSimulations);
  }
  return name;
}

std::optional<ConfigValue> parseConfigName(std::string_view name) {
  const auto tokens = split(name, '-');
  std::size_t next = 0;
  const auto peek = [&]() -> std::string_view {
    return next < tokens.size() ? tokens[next] : std::string_view{};
  };

  ConfigValue value{};

  const auto material = entryOf(MaterialNames, peek());
  if (!material.has_value()) {
    return std::nullopt;
  }
  value.materialType = material.value();
  ++next;

  const auto solver = entryOf(SolverNames, peek());
  if (!solver.has_value()) {
    return std::nullopt;
  }
  value.solver = solver.value();
  ++next;

  if (const auto mechanisms = prefixedNumber(peek(), 'm'); mechanisms.has_value()) {
    value.relaxationMechanisms = mechanisms.value();
    ++next;
  }

  const auto order = prefixedNumber(peek(), 'o');
  if (!order.has_value()) {
    return std::nullopt;
  }
  value.convergenceOrder = order.value();
  ++next;

  const auto precision = precisionOf(peek());
  if (!precision.has_value()) {
    return std::nullopt;
  }
  value.precision = precision.value();
  ++next;

  const auto drQuadRule = entryOf(DRQuadRuleNames, peek());
  if (!drQuadRule.has_value()) {
    return std::nullopt;
  }
  value.drQuadRule = drQuadRule.value();
  ++next;

  value.numSimulations = 1;
  if (const auto simulations = prefixedNumber(peek(), 's'); simulations.has_value()) {
    value.numSimulations = simulations.value();
    ++next;
  }

  if (next != tokens.size()) {
    return std::nullopt;
  }

  // every value has exactly one name; any other spelling does not name it
  if (configName(value) != name) {
    return std::nullopt;
  }
  return value;
}

std::string describeConfig(const ConfigValue& value) {
  std::ostringstream stream;
  stream << "Material: " << materialTypeName(value.materialType) << "\n";
  stream << "Solver: " << solverTypeName(value.solver) << "\n";
  stream << "Relaxation mechanisms: " << value.relaxationMechanisms << "\n";
  stream << "Convergence order: " << value.convergenceOrder << "\n";
  stream << "Precision: " << stringRealType(value.precision) << "\n";
  stream << "Dynamic rupture quadrature rule: " << drQuadRuleName(value.drQuadRule) << "\n";
  stream << "Fused simulations: " << value.numSimulations << "\n";
  return stream.str();
}

} // namespace seissol
