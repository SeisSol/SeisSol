#ifndef SEISSOL_COMMON_CELLCONFIG_HPP_
#define SEISSOL_COMMON_CELLCONFIG_HPP_

#include "Model/common_datastructures.hpp"
#include <cstddef>
#include <string>

namespace seissol {
// enum class MaterialType { Static, Acoustic, Elastic, ViscoElastic, PoroElastic, Anisotropic };

enum class PrecisionType {
  // TODO: add half precision(s) here?
  Float,
  Double
};

template <PrecisionType precision>
struct PrecisionToType {};
template <>
struct PrecisionToType<PrecisionType::Float> {
  using Type = float;
};
template <>
struct PrecisionToType<PrecisionType::Double> {
  using Type = double;
};
template <typename RealType>
struct PrecisionFromType {};
template <>
struct PrecisionFromType<float> {
  static constexpr PrecisionType Precision = PrecisionType::Float;
  static inline const std::string Text = "float";
};
template <>
struct PrecisionFromType<double> {
  static constexpr PrecisionType Precision = PrecisionType::Double;
  static inline const std::string Text = "double";
};

struct CellConfigT {
  const seissol::model::MaterialType material;
  const std::size_t mechanisms;
  const PrecisionType precision;
  const int convergenceOrder;
  const bool plasticity;

  bool operator==(const CellConfigT& other) {
    return material == other.material && mechanisms == other.mechanisms &&
           precision == other.precision && convergenceOrder == other.convergenceOrder &&
           plasticity == other.plasticity;
  }
};

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
