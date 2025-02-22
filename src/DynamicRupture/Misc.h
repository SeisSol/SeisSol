// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_MISC_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_MISC_H_

#include "Geometry/MeshDefinition.h"
#include "Kernels/Precision.h"

#include "generated_code/init.h"
#include <Initializer/Parameters/DRParameters.h>
#include <cmath>
#include <string>
#include <tuple>
#include <type_traits>

#include "Common/Marker.h"

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
static constexpr inline size_t NumPaddedPoints = leadDim<init::QInterpolated>();
static constexpr inline size_t NumQuantities = misc::dimSize<init::QInterpolated, 1>();

/**
 * Constants for Thermal Pressurization
 */
static constexpr size_t NumTpGridPoints = 60;
static constexpr double TpLogDz = 0.3;
static constexpr double TpMaxWaveNumber = 10.0;

/**
 * Number of gauss points on an element surface.
 */
static constexpr unsigned int NumBoundaryGaussPoints = init::QInterpolated::Shape[0];

template <class TupleT, class F, std::size_t... I>
constexpr F forEachImpl(TupleT&& tuple, F&& functor, std::index_sequence<I...> /*unused*/) {
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
#pragma omp declare simd
template <size_t Exp, typename T>
SEISSOL_HOSTDEVICE inline auto power(T base) -> T {
  T result = static_cast<T>(1.0);
  for (size_t i = 0; i < Exp; ++i) {
    result *= base;
  }
  return result;
}

#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> square(T t) {
  return t * t;
}

/**
 * Computes a squared sum of an N-dimensional vector
 * @return magnitude of the vector
 */
#pragma omp declare simd
template <typename T, typename... Tn>
SEISSOL_HOSTDEVICE inline T square(T t1, Tn... tn) {
  return square(t1) + square(tn...);
}

/**
 * Computes the magnitude of an N-dimensional vector
 * @return magnitude of the vector
 */
#pragma omp declare simd
template <typename T, typename... Tn>
SEISSOL_HOSTDEVICE inline T magnitude(T t1, Tn... tn) {
  return std::sqrt(square(t1) + square(tn...));
}

#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE inline T clamp(T value, T minval, T maxval) {
#ifdef __HIP__
  return std::max(minval, std::min(maxval, value));
#else
  return std::clamp(value, minval, maxval);
#endif
}

/**
 * Create strike and dip unit vectors give a fault normal vector
 * Note: equations are explained in documentation -> left-lateral-right-lateral-normal-reverse
 * @param normal
 * @param strike
 * @param dip
 */
void computeStrikeAndDipVectors(const VrtxCoords normal, VrtxCoords strike, VrtxCoords dip);

std::string frictionLawName(seissol::initializer::parameters::FrictionLawType type);

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

#endif // SEISSOL_SRC_DYNAMICRUPTURE_MISC_H_
