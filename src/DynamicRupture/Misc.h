#ifndef SEISSOL_DR_MISC_H
#define SEISSOL_DR_MISC_H

#include <cmath>
#include <generated_code/init.h>
#include <stdexcept>

#include "Kernels/precision.hpp"

namespace seissol::dr::misc {

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
constexpr int numQuantities = misc::dimSize<init::QInterpolated, 1>();
/**
 * Number of gauss points on an element surface.
 */
static constexpr unsigned int numberOfBoundaryGaussPoints =
    (CONVERGENCE_ORDER + 1) * (CONVERGENCE_ORDER + 1);

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
template <size_t exp>
double power(double base) {
  double result = 1.0;
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
real magnitude(real x, real y);

/**
 * Computes the arcus sinus hyperbolicus of x.
 * Note: precision has to be double, otherwise we would loose too much precision.
 * @param x
 * @return asinh(x)
 */
double asinh(double x);
} // namespace seissol::dr::misc

#endif // SEISSOL_DR_MISC_H
