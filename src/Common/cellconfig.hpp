#ifndef SEISSOL_COMMON_CELLCONFIG_HPP_
#define SEISSOL_COMMON_CELLCONFIG_HPP_

#include "Model/common_datastructures.hpp"
#include <cstddef>
#include <string>
#include "precision.hpp"

namespace seissol {
// TODO(David): once C++20 hits, make CellConfig implement a concept
// dummy template type
template <typename MaterialTP, typename RealTP, int ConvergenceOrderP, bool PlasticityP>
struct CellConfig {
  using MaterialT = MaterialTP;
  using RealT = RealTP;
  static constexpr bool Plasticity = PlasticityP;
  static constexpr int ConvergenceOrder = ConvergenceOrderP;

  static constexpr CellConfigT cellConfig() {
    return CellConfigT{MaterialT::Type,
                       MaterialT::Mechanisms,
                       PrecisionFromType<RealT>::Precision,
                       ConvergenceOrder,
                       Plasticity};
  }
};
} // namespace seissol

#endif
