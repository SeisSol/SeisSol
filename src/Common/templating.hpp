#pragma once

#include <type_traits>

namespace seissol {
// (the variadic pattern matching/type unpacking is partially inspired by
// https://stackoverflow.com/questions/66944744/syntax-to-unpack-tuple-on-parameter-pack-and-variadic-template
// )

// Also, say hello to a bit of templating. Just a tiny bit.

// transforms the elements of a variadic type
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

// changes the underlying variadic type
template <template <typename...> typename TargetT, typename OriginalT>
struct ChangeVariadic {};

template <template <typename...> typename TargetT,
          template <typename...>
          typename VariadicT,
          typename... Args>
struct ChangeVariadic<TargetT, VariadicT<Args...>> {
  using Result = TargetT<Args...>;
};

template <template <typename...> typename TargetT, typename OriginalT>
using ChangeVariadicT = typename ChangeVariadic<TargetT, OriginalT>::Result;

// removes type duplicates inside a variadic type
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

// a dummy type to contain variadic arguments
template <typename... Args>
struct VariadicContainer {};

// concatenates two variadic types (with the same container type)
template <typename T1, typename T2>
struct ConcatVariadic {};

template <template <typename...> typename VariadicT, typename... Args1, typename... Args2>
struct ConcatVariadic<VariadicT<Args1...>, VariadicT<Args2...>> {
  using Type = VariadicT<Args1..., Args2...>;
};

template <typename T1, typename T2>
using ConcatVariadicT = ConcatVariadic<T1, T2>;

template <typename PreT, typename T>
struct PrependVariadic {};

template <typename PreT, template <typename...> typename VariadicT, typename... Args>
struct PrependVariadic<PreT, VariadicT<Args...>> {
  using Type = VariadicT<PreT, Args...>;
};

template <typename PreT, typename T>
using PrependVariadicT = PrependVariadic<PreT, T>;

template <typename PostT, typename T>
struct AppendVariadic {};

template <typename PostT, template <typename...> typename VariadicT, typename... Args>
struct AppendVariadic<PostT, VariadicT<Args...>> {
  using Type = VariadicT<Args..., PostT>;
};

template <typename PostT, typename T>
using AppendVariadicT = AppendVariadic<PostT, T>;

} // namespace seissol
