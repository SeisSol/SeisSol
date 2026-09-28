// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONFIGDISPATCH_H_
#define SEISSOL_SRC_COMMON_CONFIGDISPATCH_H_

#include "Common/ConfigRegistry.h"
#include "Config.h"

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <variant>

namespace seissol {

/**
 * @brief The configurations built into the executable as types, in the order of their ids.
 *
 * Code that depends on a configuration as a type gets there from a `ConfigId` through
 * `dispatchConfig`, and back through `configIdOf`.
 */
using ConfigVariant = std::variant<Config>;

namespace internal {

template <typename Cfg, std::size_t Index = 0>
constexpr ConfigId configIdOf() {
  static_assert(Index < std::variant_size_v<ConfigVariant>,
                "The configuration is not built into the executable.");
  if constexpr (std::is_same_v<std::variant_alternative_t<Index, ConfigVariant>, Cfg>) {
    return Index;
  } else {
    return configIdOf<Cfg, Index + 1>();
  }
}

template <std::size_t... Indices>
constexpr std::array<ConfigVariant, sizeof...(Indices)>
    configVariants(std::index_sequence<Indices...> /*indices*/) {
  return {ConfigVariant(std::in_place_index<Indices>)...};
}

} // namespace internal

/// The id of a configuration type.
template <typename Cfg>
constexpr ConfigId configIdOf() {
  return internal::configIdOf<Cfg>();
}

/**
 * @brief Calls `function` with the configuration type that `id` stands for.
 *
 * The configuration is passed as a value-initialized object of its type; `decltype` recovers the
 * type from it.
 */
template <typename F>
decltype(auto) dispatchConfig(ConfigId id, F&& function) {
  static constexpr auto Variants =
      internal::configVariants(std::make_index_sequence<std::variant_size_v<ConfigVariant>>());
  return std::visit(std::forward<F>(function), Variants.at(id));
}

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONFIGDISPATCH_H_
