#ifndef SEISSOL_DR_MISC_H
#define SEISSOL_DR_MISC_H

#include "Geometry/MeshDefinition.h"
#include "Kernels/precision.hpp"

#include <cmath>
#include <generated_code/init.h>
#include <stdexcept>
#include <tuple>

namespace seissol::dr::misc {
// TODO: this can be moved to yateto headers
template <typename Tensor, int Dim>
constexpr size_t dimSize() noexcept {
  return Tensor::Stop[Dim] - Tensor::Start[Dim];
}

template <typename Tensor>
constexpr size_t leadDim() noexcept {
  return dimSize<Tensor, 0>();
}

/**
 * Number of gauss points padded to match the vector register length.
 */
static constexpr inline size_t numPaddedPoints = leadDim<init::QInterpolated>();
static constexpr inline size_t numQuantities = misc::dimSize<init::QInterpolated, 1>();

/**
 * Constants for Thermal Pressurization
 */
static constexpr size_t numberOfTPGridPoints = 60;
static constexpr real tpLogDz = 0.3;
static constexpr real tpMaxWaveNumber = 10.0;

/**
 * Number of gauss points on an element surface.
 */
static constexpr unsigned int numberOfBoundaryGaussPoints = init::QInterpolated::Shape[0];

template <class TupleT, class F, std::size_t... I>
constexpr F forEachImpl(TupleT&& tuple, F&& functor, std::index_sequence<I...>) {
  return (void)std::initializer_list<int>{
             (std::forward<F>(functor)(std::get<I>(std::forward<TupleT>(tuple)), I), 0)...},
         functor;
}

template <typename TupleT, typename F>
constexpr F forEach(TupleT&& tuple, F&& functor) {
  return forEachImpl(
      std::forward<TupleT>(tuple),
      std::forward<F>(functor),
      std::make_index_sequence<std::tuple_size<std::remove_reference_t<TupleT>>::value>{});
}
/**
 * Compute base^exp
 * Note: precision has to be double, otherwise we would loose too much precision.
 * @param base
 * @return
 */
template <size_t exp, typename T>
inline auto power(T base) -> T {
  T result = static_cast<T>(1.0);
  for (size_t i = 0; i < exp; ++i) {
    result *= base;
  }
  return result;
}

/**
 * Computes the magnitude of the vector (x, y)
 * @param x First component of the vector
 * @param y Second component of the vector
 * @return magnitude of the vector
 */
inline real magnitude(real x, real y) { return std::sqrt(x * x + y * y); }

/**
 * Computes the magnitude of the vector (x, y)
 * @param x First component of the vector
 * @param y Second component of the vector
 * @return magnitude of the vector
 */
inline real magnitude(real x, real y, real z) { return std::sqrt(x * x + y * y + z * z); }

/**
 * Computes the arcus sinus hyperbolicus of x.
 * Note: precision has to be double, otherwise we would loose too much precision.
 * @param x
 * @return asinh(x)
 */
inline double asinh(double x) { return std::log(x + std::sqrt(x * x + 1.0)); }

/**
 * Create strike and dip unit vectors give a fault normal vector
 * Note: equations are explained in documentation -> left-lateral-right-lateral-normal-reverse
 * @param normal
 * @param strike
 * @param dip
 */
void computeStrikeAndDipVectors(const VrtxCoords normal, VrtxCoords strike, VrtxCoords dip);

namespace quantity_indices {
/**
 * Defines the indices under which one can find a specific quantity.
 * U, V, W: Velocities in x, y, z direction.
 * N, T1, T2: traction in normal and fault aligned directions.
 * XX, YY, ZZ, XY, YZ, XZ: Stress in cartesian coordinates
 * Use as:
 * ```
 * using namepace dr::misc::quantity_indices;
 * real quantities[9];
 * real normalStress = quantities[N];
 * ```
 * */
enum QuantityIndices : size_t {
  U = 6,
  V = 7,
  W = 8,
  N = 0,
  T1 = 3,
  T2 = 5,
  XX = 0,
  YY = 1,
  ZZ = 2,
  XY = 3,
  YZ = 4,
  XZ = 5,
};
} // namespace quantity_indices
} // namespace seissol::dr::misc

#endif // SEISSOL_DR_MISC_H
