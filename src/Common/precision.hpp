#pragma once

#include <string>

namespace seissol {

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

} // namespace seissol
