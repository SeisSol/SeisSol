#pragma once

#include <type_traits>

namespace seissol {
// (partially inspired by
// https://stackoverflow.com/questions/66944744/syntax-to-unpack-tuple-on-parameter-pack-and-variadic-template
// )

// Also, say hello to a bit of templating. Just a tiny bit.

template <template <typename> typename ElementTransform, typename OriginalT>
struct TransformVariadic {};

template <template <typename> typename ElementTransform,
          template <typename...>
          typename VariadicT,
          typename... Args>
struct TransformVariadic<ElementTransform, VariadicT<Args...>> {
  using Result = VariadicT<ElementTransform<Args>...>;
};

template <template <typename> typename ElementTransform, typename OriginalT>
using TransformVariadicT = typename TransformVariadic<ElementTransform, OriginalT>::Result;

template <typename OriginalT>
struct RemoveDuplicateVariadic {};

template <template <typename...> typename VariadicT, typename... Args>
struct RemoveDuplicateVariadic<VariadicT<Args...>> {
  template <typename Head, typename... Rest>
  constexpr static bool containsHead() {
    return (false || ... || std::is_same_v<Head, Rest>);
  }
  template <typename Head, typename T>
  struct VariadicPrepend {};
  template <typename Head, typename... Rest>
  struct VariadicPrepend<Head, VariadicT<Rest...>> {
    using Result = VariadicT<Head, Rest...>;
  };
  template <typename Head, typename... Rest>
  struct Intermediate {
    using PreResult = typename Intermediate<Rest...>::Result;
    using Result = std::conditional_t<containsHead<Head, Rest...>(),
                                      PreResult,
                                      typename VariadicPrepend<Head, PreResult>::Result>;
  };
  template <typename Head>
  struct Intermediate<Head> {
    using Result = VariadicT<Head>;
  };

  using Result = typename Intermediate<Args...>::Result;
};

template <typename OriginalT>
using RemoveDuplicateVariadicT = typename RemoveDuplicateVariadic<OriginalT>::Result;

} // namespace seissol
