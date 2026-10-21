// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_TESTCONFIGS_H_
#define SEISSOL_TESTS_TESTCONFIGS_H_

// The configurations built into the executable, for the tests that run once per configuration
// (TEST_CASE_TEMPLATE with SEISSOL_CONFIG_TYPES, or TEST_CASE_TEMPLATE_APPLY with one of the lists
// below). Include it before the first test of a configuration, so that doctest lists each instance
// with the canonical name of its configuration.

#include <doctest.h>

#include "Common/ConfigValue.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Model/MaterialType.h"

#include <string>
#include <tuple>
#include <type_traits>

namespace doctest::detail {
#define SEISSOL_TEST_CONFIG_NAME(Cfg)                                                              \
  template <>                                                                                      \
  inline const char* type_to_string<Cfg>() {                                                       \
    static const std::string Name = "<" + ::seissol::configName(Cfg::Value) + ">";                 \
    return Name.c_str();                                                                           \
  }
SEISSOL_FOR_EACH_CONFIG(SEISSOL_TEST_CONFIG_NAME)
#undef SEISSOL_TEST_CONFIG_NAME
} // namespace doctest::detail

namespace seissol::unit_test {

namespace internal {

template <typename Tuple>
struct TupleTail;

template <typename Head, typename... Tail>
struct TupleTail<std::tuple<Head, Tail...>> {
  using Type = std::tuple<Tail...>;
};

template <typename Cfg, typename Tuple>
struct TuplePrepend;

template <typename Cfg, typename... Cfgs>
struct TuplePrepend<Cfg, std::tuple<Cfgs...>> {
  using Type = std::tuple<Cfg, Cfgs...>;
};

template <template <typename> typename Pred, typename... Cfgs>
struct ConfigFilter {
  using Type = std::tuple<>;
};

template <template <typename> typename Pred, typename Cfg, typename... Rest>
struct ConfigFilter<Pred, Cfg, Rest...> {
  using RestType = typename ConfigFilter<Pred, Rest...>::Type;
  using Type =
      std::conditional_t<Pred<Cfg>::value, typename TuplePrepend<Cfg, RestType>::Type, RestType>;
};

template <model::MaterialType Type>
struct MaterialIs {
  template <typename Cfg>
  using Pred = std::bool_constant<Cfg::MaterialType == Type>;
};

template <SolverType Type>
struct SolverIs {
  template <typename Cfg>
  using Pred = std::bool_constant<Cfg::Solver == Type>;
};

template <typename Cfg>
struct SupportsDynamicRupture : std::bool_constant<model::MaterialOf<Cfg>::SupportsDR> {};

#define SEISSOL_TEST_CONFIG_ENTRY(Cfg) , Cfg
using MaterialConfigs =
    TupleTail<std::tuple<void SEISSOL_FOR_EACH_MATERIAL(SEISSOL_TEST_CONFIG_ENTRY)>>::Type;
#undef SEISSOL_TEST_CONFIG_ENTRY

} // namespace internal

/// The configurations built into the executable for which `Pred<Cfg>::value` holds, as a
/// `std::tuple` for TEST_CASE_TEMPLATE_APPLY. A test applied to none of them does not exist.
template <template <typename> typename Pred>
using ConfigsWhere = typename internal::ConfigFilter<Pred, SEISSOL_CONFIG_TYPES>::Type;

/// The configurations built into the executable whose material is of the type `Type`.
template <model::MaterialType Type>
using ConfigsOfMaterial = ConfigsWhere<internal::MaterialIs<Type>::template Pred>;

/// The configurations built into the executable that advance their cells with the solver `Type`.
template <SolverType Type>
using ConfigsOfSolver = ConfigsWhere<internal::SolverIs<Type>::template Pred>;

/// The configurations built into the executable whose material supports dynamic rupture.
using DynamicRuptureConfigs = ConfigsWhere<internal::SupportsDynamicRupture>;

/// The first configuration of each material built into the executable (SEISSOL_FOR_EACH_MATERIAL);
/// for the tests of a material, `model::MaterialOf<Cfg>`.
using MaterialConfigs = internal::MaterialConfigs;

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_TESTCONFIGS_H_
