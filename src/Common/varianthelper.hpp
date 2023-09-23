#pragma once

#include <Common/templating.hpp>
#include <functional>
#include <iterator>
#include <tuple>
#include <variant>

namespace seissol {
template <typename T>
using VariantStorage = ChangeVariadicT<std::tuple, T>;

template <std::size_t Idx, typename T, typename VariantT>
constexpr std::size_t variantTypeIdInner() {
  static_assert(Idx < std::variant_size_v<VariantT>, "Type not found in variant");
  if constexpr (std::is_same_v<std::variant_alternative_t<Idx, VariantT>, T>) {
    return Idx;
  }
  variantTypeIdInner<Idx + 1, T, VariantT>();
}

template <typename T, typename VariantT>
constexpr std::size_t variantTypeId() {
  return variantTypeIdInner<0, T, VariantT>();
}

template <std::size_t Idx, typename T, typename F>
void variantContainerForAllInner(T& data, F&& visitor) {
  if constexpr (Idx < std::tuple_size_v<T>) {
    std::invoke(visitor, data);
    variantContainerForAllInner<Idx + 1>(data, std::forward<F>(visitor));
  }
}

template <typename T, typename F>
void variantContainerForAll(T& data, F&& visitor) {
  variantContainerForAllInner<0>(data, std::forward<F>(visitor));
}
} // namespace seissol
