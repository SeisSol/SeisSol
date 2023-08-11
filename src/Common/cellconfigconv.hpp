#pragma once

#include "cellconfig.hpp"
#include "configs.hpp"
#include "utils/logger.h"
#include <type_traits>
#include <variant>

namespace seissol {
template <std::size_t I>
void printSupportedConfigsSub() {
  if constexpr (I < std::variant_size_v<SupportedConfigs>) {
    using ConfigI = std::variant_alternative_t<I, SupportedConfigs>;
    logInfo() << I << ":::" << ConfigI::MaterialT::Text
              << PrecisionFromType<typename ConfigI::RealT>::Text << ConfigI::ConvergenceOrder
              << ConfigI::Plasticity;
    printSupportedConfigsSub<I + 1>();
  }
}

void printSupportedConfigs() {
  logInfo() << "The following cell configurations are supported in this build of SeisSol:";
  printSupportedConfigsSub<0>();
  logInfo() << "The end.";
}

// partially inspired by
// https://stackoverflow.com/questions/66944744/syntax-to-unpack-tuple-on-parameter-pack-and-variadic-template

// Also, say hello to a bit of templating. Just a tiny bit.

template <typename OriginalT, template <typename> typename ElementTransform>
struct TransformVariadic {};

template <template <typename> typename ElementTransform,
          template <typename...>
          typename VariadicT,
          typename... Args>
struct TransformVariadic<VariadicT<Args...>, ElementTransform> {
  using Result = VariadicT<ElementTransform<Args>...>;
};

template <typename OriginalT>
struct RemoveDuplicateVariadic {};

template <template <typename...> typename VariadicT, typename... Args>
struct RemoveDuplicateVariadic<VariadicT<Args...>> {
  template <typename Head>
  constexpr static bool containsHead() {
    return false;
  }
  template <typename Head, typename Head2, typename... Rest>
  constexpr static bool containsHead() {
    return containsHead<Head, Rest...>() || std::is_same_v<Head, Head2>;
  }
  template <typename T>
  struct VariadicPrepend {};
  template <typename... Rest>
  struct VariadicPrepend<VariadicT<Rest...>> {
    template <typename Head>
    using Prepend = VariadicT<Head, Rest...>;
  };
  template <typename Head, typename... Rest>
  struct Intermediate {
    using PreResult = typename Intermediate<Rest...>::Result;
    using Result = std::conditional_t<containsHead<Head, Rest...>(),
                                      PreResult,
                                      typename VariadicPrepend<PreResult>::template Prepend<Head>>;
  };
  template <typename Head>
  struct Intermediate<Head> {
    using Result = VariadicT<Head>;
  };

  using Result = typename Intermediate<Args...>::Result;
};

template <typename Config>
using SelectMaterial = typename Config::MaterialT;
template <typename Config>
using SelectReal = typename Config::RealT;

using SupportedMaterials =
    RemoveDuplicateVariadic<TransformVariadic<SupportedConfigs, SelectMaterial>::Result>::Result;
using SupportedReals =
    RemoveDuplicateVariadic<TransformVariadic<SupportedConfigs, SelectReal>::Result>::Result;

constexpr SupportedConfigs defaultConfig(bool plasticity) {
  if (plasticity) {
    return SupportedConfigs(CellConfig<seissol::model::Material_t, real, ConvergenceOrder, true>());
  } else {
    return SupportedConfigs(
        CellConfig<seissol::model::Material_t, real, ConvergenceOrder, false>());
  }
}

template <typename OriginalT>
struct DeclareVariadic {};

template <template <typename...> typename VariadicT, typename... Args>
struct DeclareVariadic<VariadicT<Args...>> {
  void dummy(Args... args) {}
};

template <template <typename> typename Struct>
using DeclareForAllConfigs = DeclareVariadic<TransformVariadic<SupportedConfigs, Struct>>;

} // namespace seissol
