#pragma once

#include "cellconfigconv.hpp"
#include <Common/configtensor.hpp>
#include <Kernels/common.hpp>
#include <Model/common_datastructures.hpp>
#include <Model/plasticity.hpp>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include "utils/logger.h"

namespace {
template <std::size_t I>
static constexpr std::array<seissol::SupportedConfigs,
                            std::variant_size_v<seissol::SupportedConfigs>>
    configArray(std::array<seissol::SupportedConfigs,
                           std::variant_size_v<seissol::SupportedConfigs>>&& array) {
  if constexpr (I < std::variant_size_v<seissol::SupportedConfigs>) {
    array[I] = std::variant_alternative_t<I, seissol::SupportedConfigs>();
    return configArray<I + 1>(
        std::forward<
            std::array<seissol::SupportedConfigs, std::variant_size_v<seissol::SupportedConfigs>>>(
            array));
  } else {
    return array;
  }
}

template <std::size_t I>
void printSupportedConfigsSub() {
  if constexpr (I < std::variant_size_v<seissol::SupportedConfigs>) {
    using ConfigI = std::variant_alternative_t<I, seissol::SupportedConfigs>;
    logInfo() << I << ":::" << ConfigI::MaterialT::Text
              << seissol::PrecisionFromType<typename ConfigI::RealT>::Text
              << ConfigI::ConvergenceOrder << ConfigI::Plasticity;
    printSupportedConfigsSub<I + 1>();
  }
}
} // namespace

namespace seissol {
template <typename Config>
struct VerifyTensorSizes {
  using MaterialT = typename Config::MaterialT;
  using RealT = typename Config::RealT;

  using Constants = ConfigConstants<Config>;

  static_assert(std::is_floating_point_v<RealT>,
                "CellConfig correctness check failed: RealT needs to be a floating point type.");
  static_assert(
      std::is_base_of_v<seissol::model::Material, MaterialT>,
      "CellConfig correctness check failed: MaterialT needs to subclass the Material class.");

  static_assert(std::is_same_v<RealT, real>,
                "CellConfig correctness check failed: RealT does not match the internal real type "
                "[KERNELS MISSING].");
  static_assert(Constants::DofsElaSize == Yateto<Config>::Tensor::Q::size(),
                "CellConfig correctness check failed: tensor sizes do no match [KERNELS MISSING].");
  static_assert(
      Constants::DofsAneSize ==
          Yateto<Config>::Kernel::template size<Yateto<Config>::Tensor::Qane>(),
      "CellConfig correctness check failed: anelastic tensor sizes do no match [KERNELS MISSING].");

  static_assert(::seissol::GivenNumberOfQuantities == Config::MaterialT::NumberOfQuantities,
                "CellConfig correctness check failed: The provided number of quantities (i.e. for "
                "which code has been generated) does not match the number of quantities specified "
                "by the material [KERNELS MISSING].");

  static_assert(
      Config::ConvergenceOrder <= yateto::numFamilyMembers<Yateto<Config>::Tensor::dQ>(),
      "CellConfig correctness check failed: Too high convergence order for generated kernels.");
};

void printSupportedConfigs() {
  logInfo() << "The following cell configurations are supported in this build of SeisSol:";
  printSupportedConfigsSub<0>();
  logInfo() << "The end.";
}

const std::array<SupportedConfigs, std::variant_size_v<SupportedConfigs>> ConfigInstances =
    configArray<0>({});
} // namespace seissol

namespace seissol::_definitions {
const seissol::DeclareForAllConfigs<seissol::VerifyTensorSizes> declVerify;
} // namespace seissol::_definitions
