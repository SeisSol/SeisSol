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
    logInfo() << I << ":::" << ConfigI::MaterialT::Text << PrecisionFromType<typename ConfigI::RealT>::Text << ConfigI::ConvergenceOrder << ConfigI::Plasticity;
    printSupportedConfigsSub<I+1>();
  }
}

void printSupportedConfigs() {
    logInfo() << "The following cell configurations are supported in this build of SeisSol:";
    printSupportedConfigsSub<0>();
    logInfo() << "The end.";
}

template <std::size_t I>
SupportedConfigs fromStructSub(const CellConfigT& config) {
  if constexpr (I < std::variant_size_v<SupportedConfigs>) {
    using ConfigI = std::variant_alternative_t<I, SupportedConfigs>;
    if (ConfigI::cellConfig() == config) {
      return SupportedConfigs(ConfigI());
    } else {
      return fromStructSub<I + 1>(config);
    }
  } else {
    // TODO(David): make more descriptive
    logError() << "Unknown cell configuration.";
    throw std::runtime_error("Unknown cell configuration.");
  }
}

SupportedConfigs configFromStruct(const CellConfigT& config) { return fromStructSub<0>(config); }

constexpr CellConfigT configToStruct(const SupportedConfigs& config) {
  return std::visit(
      [&](auto&& elem) {
        using ConfigT = std::decay_t<decltype(elem)>;
        return ConfigT::cellConfig();
      },
      config);
}

// partially inspired by https://stackoverflow.com/questions/66944744/syntax-to-unpack-tuple-on-parameter-pack-and-variadic-template

template<typename OriginalT, template<typename> typename ElementTransform>
struct TransformVariadic {
};

template<template<typename> typename ElementTransform, template<typename...> typename VariadicT, typename ...Args>
struct TransformVariadic<VariadicT<Args...>, ElementTransform> {
  using Result = VariadicT<ElementTransform<Args>...>;
};

template<typename Config>
using SelectMaterial = typename Config::MaterialT;
template<typename Config>
using SelectReal = typename Config::RealT;

using SupportedMaterials = TransformVariadic<SupportedConfigs, SelectMaterial>::Result;
using SupportedReals = TransformVariadic<SupportedConfigs, SelectReal>::Result;

constexpr SupportedConfigs defaultConfig(bool plasticity) {
  if (plasticity) {
    return SupportedConfigs(CellConfig<seissol::model::Material_t, real, ConvergenceOrder, true>());
  } else {
    return SupportedConfigs(
        CellConfig<seissol::model::Material_t, real, ConvergenceOrder, false>());
  }
}

template<typename OriginalT>
struct DeclareVariadic {
};

template<template<typename...> typename VariadicT, typename ...Args>
struct DeclareVariadic<VariadicT<Args...>> {
  void dummy(Args... args) {}
};

template<template<typename> typename Struct>
using DeclareForAllConfigs = DeclareVariadic<TransformVariadic<SupportedConfigs, Struct>>;

} // namespace seissol
