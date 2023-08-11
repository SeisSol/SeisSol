#pragma once

#include <Kernels/common.hpp>
#include <Model/plasticity.hpp>
#include <cstddef>
#include "constants.hpp"

namespace seissol {
template <typename Config>
struct ConfigConstants {
  ConfigConstants() = delete; // static class

  using RealT = typename Config::RealT;
  using MaterialT = typename Config::MaterialT;

  static constexpr std::size_t ElaSize = MaterialT::NumberOfQuantities;
  static constexpr std::size_t AneSize = MaterialT::NumberPerMechanism * MaterialT::Mechanisms;

  static constexpr std::size_t DofsElaSize =
      ElaSize * seissol::kernels::NumberOfBasisFunctions(Config::ConvergenceOrder);
  static constexpr std::size_t DofsAneSize =
      AneSize * seissol::kernels::NumberOfBasisFunctions(Config::ConvergenceOrder);
  static constexpr std::size_t PStrainSize =
      seissol::model::PlasticityData<RealT>::NumberOfQuantities *
      seissol::kernels::NumberOfAlignedBasisFunctions(Config::ConvergenceOrder);

  static constexpr std::size_t TensorSizeQ = DofsElaSize * Config::ConvergenceOrder;
  static constexpr std::size_t TensorSizeI = DofsElaSize;
  static constexpr std::size_t TensorSizeQAne =
      ZeroLengthArrayHandler(DofsAneSize * Config::ConvergenceOrder);
  static constexpr std::size_t TensorSizeIAne = ZeroLengthArrayHandler(DofsAneSize);
};
} // namespace seissol
