#pragma once

namespace seissol {

enum class RealType { F32, F64 };

template <RealType P>
struct RealTypeWrapper {};

template <>
struct RealTypeWrapper<RealType::F32> {
  using Type = float;
};

template <>
struct RealTypeWrapper<RealType::F64> {
  using Type = double;
};

template <RealType P>
using RealT = typename RealTypeWrapper<P>::Type;

} // namespace seissol
