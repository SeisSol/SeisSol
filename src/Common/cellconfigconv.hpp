#pragma once

#include "templating.hpp"
#include "cellconfig.hpp"
#include "configs.hpp"
#include <Kernels/common.hpp>
#include "Model/plasticity.hpp"
#include <type_traits>
#include <variant>

namespace seissol {
void printSupportedConfigs();

template <typename Config>
using SelectMaterial = typename Config::MaterialT;
template <typename Config>
using SelectReal = typename Config::RealT;

using SupportedMaterials =
    RemoveDuplicateVariadicT<TransformVariadicT<SelectMaterial, SupportedConfigs>>;
using SupportedReals = RemoveDuplicateVariadicT<TransformVariadicT<SelectReal, SupportedConfigs>>;

template <typename OriginalT>
struct DeclareVariadic {};

template <template <typename...> typename VariadicT, typename... Args>
struct DeclareVariadic<VariadicT<Args...>> {
  void dummy(Args... args) {}
};

template <template <typename> typename Struct>
using DeclareForAllConfigs = DeclareVariadic<TransformVariadicT<Struct, SupportedConfigs>>;

extern const std::array<SupportedConfigs, std::variant_size_v<SupportedConfigs>> ConfigInstances;

template <typename Config>
struct ConfigConstants {
  ConfigConstants() = delete; // static class

  using RealT = typename Config::RealT;
  using MaterialT = typename Config::MaterialT;

  static constexpr std::size_t ElaSize = MaterialT::NumberOfQuantities;
  static constexpr std::size_t AneSize = MaterialT::NumberPerMechanism * MaterialT::Mechanisms;

  static constexpr std::size_t DofsElaSize =
      ElaSize * seissol::kernels::NumberOfBasisFunctions(Config::ConvergenceOrder);
  static constexpr std::size_t DerivativesElaSize = DofsElaSize * Config::ConvergenceOrder;
  static constexpr std::size_t DofsAneSize =
      AneSize * seissol::kernels::NumberOfBasisFunctions(Config::ConvergenceOrder);
  static constexpr std::size_t PStrainSize =
      seissol::model::PlasticityData<RealT>::NumberOfQuantities *
      seissol::kernels::NumberOfAlignedBasisFunctions<RealT>(Config::ConvergenceOrder);

  static constexpr std::size_t TensorSizeQ = DofsElaSize * Config::ConvergenceOrder;
  static constexpr std::size_t TensorSizeI = DofsElaSize;
  static constexpr std::size_t TensorSizeQAne =
      ZeroLengthArrayHandler(DofsAneSize * Config::ConvergenceOrder);
  static constexpr std::size_t TensorSizeIAne = ZeroLengthArrayHandler(DofsAneSize);

  // TODO(David): maybe remove those again
  static constexpr std::size_t TensorSizew = MaterialT::Mechanisms;
  static constexpr std::size_t TensorSizeW = MaterialT::Mechanisms * MaterialT::Mechanisms;
  static constexpr std::size_t TensorSizeE = ElaSize * AneSize;
  static constexpr std::size_t TensorSizeET =
      MaterialT::NumberOfQuantities * MaterialT::NumberOfQuantities;
};

} // namespace seissol
