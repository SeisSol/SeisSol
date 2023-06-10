#pragma once

#include "cellconfig.hpp"
#include "configs.hpp"
#include "utils/logger.h"
#include <typeinfo>

namespace seissol {
template <std::size_t I>
void printSupportedConfigsSub() {
  if constexpr (I < std::variant_size_v <SupportedConfigs>) {
    using ConfigI = std::variant_alternative_t<I, SupportedConfigs>;
    logInfo() << I << ":::" << typeid(ConfigI).name();
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

constexpr SupportedConfigs defaultConfig(bool plasticity) {
  if (plasticity) {
    return SupportedConfigs(CellConfig<seissol::model::Material_t, real, ConvergenceOrder, true>());
  } else {
    return SupportedConfigs(
        CellConfig<seissol::model::Material_t, real, ConvergenceOrder, false>());
  }
}
} // namespace seissol
