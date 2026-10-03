// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONFIGVALUE_H_
#define SEISSOL_SRC_COMMON_CONFIGVALUE_H_

#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Model/MaterialType.h"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace seissol {

/**
 * @brief The choices that select the kernels of a build, as a runtime value.
 *
 * A compile-time configuration (the generated `Config`) carries one of these as its `Value` and
 * takes each of its constants from the corresponding field, so both describe the same thing. A
 * value can be stored, compared and looked up, which is what code needs that does not know the
 * configuration as a type.
 *
 * What the choices imply for the generated tensors (their sizes, the number of points on a face,
 * ...) is not part of the value: that depends on the target as well.
 *
 * A value-initialized `ConfigValue` names no configuration.
 */
struct ConfigValue {
  std::size_t convergenceOrder{};
  std::size_t relaxationMechanisms{};
  model::MaterialType materialType{};
  RealType precision{};
  SolverType solver{};
  DRQuadRuleType drQuadRule{};
  std::size_t numSimulations{};
};

constexpr bool operator==(const ConfigValue& lhs, const ConfigValue& rhs) {
  return lhs.convergenceOrder == rhs.convergenceOrder &&
         lhs.relaxationMechanisms == rhs.relaxationMechanisms &&
         lhs.materialType == rhs.materialType && lhs.precision == rhs.precision &&
         lhs.solver == rhs.solver && lhs.drQuadRule == rhs.drQuadRule &&
         lhs.numSimulations == rhs.numSimulations;
}

constexpr bool operator!=(const ConfigValue& lhs, const ConfigValue& rhs) { return !(lhs == rhs); }

/// The name of a material in a configuration name, e.g. `elastic`.
std::string_view materialTypeName(model::MaterialType type);

/// The name of a solver in a configuration name, e.g. `linearck`.
std::string_view solverTypeName(SolverType type);

/// The name of a quadrature rule for dynamic rupture faces in a configuration name, e.g. `stroud`.
std::string_view drQuadRuleName(DRQuadRuleType rule);

/**
 * @brief The canonical name of a configuration.
 *
 * It lists every field of the value in a fixed order, as
 *
 *     <material>-<solver>[-m<mechanisms>]-o<order>-<precision>-<quadrature rule>[-s<simulations>]
 *
 * e.g. `elastic-linearck-o6-f64-stroud`. The mechanisms are left out when there are none, the
 * simulations when there is only one. Distinct values have distinct names, so a name can stand
 * for a configuration outside the executable, e.g. in a file.
 */
std::string configName(const ConfigValue& value);

/**
 * @brief The value a canonical name stands for.
 *
 * Returns nothing for every string that `configName` does not produce, including other spellings
 * of a valid configuration (such as `o06` or an explicit `s1`).
 */
std::optional<ConfigValue> parseConfigName(std::string_view name);

/// A description of the configuration for humans; one field per line.
std::string describeConfig(const ConfigValue& value);

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONFIGVALUE_H_
