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

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
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
using ConfigVariant = std::variant<SEISSOL_CONFIG_TYPES>;

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

template <template <typename> typename T, typename Variant>
struct PerConfigValues;

template <template <typename> typename T, typename... Cfgs>
struct PerConfigValues<T, std::variant<Cfgs...>> {
  using Type = std::tuple<T<Cfgs>...>;
};

template <template <typename> typename T, typename Variant>
struct OfConfigVariant;

template <template <typename> typename T, typename... Cfgs>
struct OfConfigVariant<T, std::variant<Cfgs...>> {
  using Type = std::variant<T<Cfgs>...>;
};

constexpr bool forEachConfigListsVariant() {
  std::size_t index = 0;
  bool listed = true;
#define SEISSOL_CONFIG_LISTED(Cfg) listed = listed && configIdOf<Cfg>() == index++;
  SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_LISTED)
#undef SEISSOL_CONFIG_LISTED
  return listed && index == std::variant_size_v<ConfigVariant>;
}

} // namespace internal

static_assert(internal::forEachConfigListsVariant(),
              "SEISSOL_FOR_EACH_CONFIG lists the configurations of ConfigVariant, in its order.");

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

/// Calls `function` as `dispatchConfig` does, once for every configuration built into the
/// executable, in the order of their ids.
template <typename F>
void forEachConfig(const F& function) {
  for (ConfigId id = 0; id < std::variant_size_v<ConfigVariant>; ++id) {
    dispatchConfig(id, function);
  }
}

/// The largest convergence order among the configurations built into the executable.
constexpr std::size_t maxBuiltConvergenceOrder() {
  std::size_t order = 0;
#define SEISSOL_CONFIG_ORDER(Cfg) order = std::max(order, Cfg::ConvergenceOrder);
  SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_ORDER)
#undef SEISSOL_CONFIG_ORDER
  return order;
}

/**
 * @brief One `T<Cfg>` for every configuration `Cfg` built into the executable.
 *
 * Holds what exists once per configuration, such as the global data of its kernels.
 */
template <template <typename> typename T>
class PerConfig {
  public:
  template <typename Cfg>
  [[nodiscard]] T<Cfg>& get() {
    return std::get<configIdOf<Cfg>()>(values_);
  }

  template <typename Cfg>
  [[nodiscard]] const T<Cfg>& get() const {
    return std::get<configIdOf<Cfg>()>(values_);
  }

  private:
  typename internal::PerConfigValues<T, ConfigVariant>::Type values_;
};

/**
 * @brief A `T<Cfg>` for one of the configurations `Cfg` built into the executable, at the index of
 * its id.
 *
 * Holds what belongs to an object of some configuration, such as the transformations of a fault
 * face.
 */
template <template <typename> typename T>
using ConfigVariantOf = typename internal::OfConfigVariant<T, ConfigVariant>::Type;

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONFIGDISPATCH_H_
