#ifndef SEISSOL_COMMON_CELLCONFIG_HPP_
#define SEISSOL_COMMON_CELLCONFIG_HPP_

#include <cstddef>
#include "Model/common_datastructures.hpp"

// shortening of arguments... Not ideal. But it works
#define CELLCONFIG_TEMPLATE                                                                        \
  typename MaterialT, typename RealT, int ConvergenceOrder, bool Plasticity
#define CELLCONFIG_TEMPLATE_FWD MaterialT, RealT, ConvergenceOrder, Plasticity

#define CELLCONFIG_TEMPLATE1 typename RealT, int ConvergenceOrder, bool Plasticity

namespace seissol {
// enum class MaterialType { Static, Acoustic, Elastic, ViscoElastic, PoroElastic, Anisotropic };

enum class PrecisionType {
  // TODO: add half precision(s) here?
  F32,
  F64
};

template <PrecisionType precision>
struct PrecisionToType {};
template <>
struct PrecisionToType<PrecisionType::F32> {
  using Type = float;
};
template <>
struct PrecisionToType<PrecisionType::F64> {
  using Type = double;
};
template <typename RealType>
struct PrecisionFromType {};
template <>
struct PrecisionFromType<float> {
  static constexpr PrecisionType Precision = PrecisionType::F32;
};
template <>
struct PrecisionFromType<double> {
  static constexpr PrecisionType Precision = PrecisionType::F64;
};

struct CellConfigT {
  const seissol::model::MaterialType material;
  const std::size_t mechanisms;
  const PrecisionType real;
  const int convergenceOrder;
  const bool plasticity;

  bool operator==(const CellConfigT& other) {
    return material == other.material && mechanisms == other.mechanisms && real == other.real &&
           convergenceOrder == other.convergenceOrder && plasticity == other.plasticity;
  }
};

// TODO(David): once C++20 hits, make CellConfig implement a concept
// dummy template type
template <typename MaterialT, typename RealT, int ConvergenceOrderP, bool PlasticityP>
struct CellConfig {
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
