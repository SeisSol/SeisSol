#pragma once

#include "cellconfig.hpp"
#include "configs.hpp"
#include "utils/logger.h"
#include <type_traits>
#include <variant>
#include "templating.hpp"

namespace seissol {
void printSupportedConfigs();

template <typename Config>
using SelectMaterial = typename Config::MaterialT;
template <typename Config>
using SelectReal = typename Config::RealT;

using SupportedMaterials =
    RemoveDuplicateVariadicT<TransformVariadicT<SelectMaterial, SupportedConfigs>>;
using SupportedReals = RemoveDuplicateVariadicT<TransformVariadicT<SelectReal, SupportedConfigs>>;

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
using DeclareForAllConfigs = DeclareVariadic<TransformVariadicT<Struct, SupportedConfigs>>;

extern const std::array<SupportedConfigs, std::variant_size_v<SupportedConfigs>> ConfigInstances;

} // namespace seissol
