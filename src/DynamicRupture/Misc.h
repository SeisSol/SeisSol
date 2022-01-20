#ifndef SEISSOL_DR_MISC_H
#define SEISSOL_DR_MISC_H

#include <generated_code/init.h>
#include <stdexcept>
#include <cmath>

#include "Kernels/precision.hpp"

namespace seissol::dr::misc {
template <typename Tensor>
constexpr size_t leadDim() noexcept {
  return Tensor::Stop[0] - Tensor::Start[0];
}
static constexpr inline size_t AlignedNumGaussPoints = leadDim<init::QInterpolated>();

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

template <size_t exp>
real power(real base) {
  real result = 1.0;
  for (size_t i = 0; i < exp; ++i) {
    result *= base;
  }
  return result;
}

real magnitude(real x, real y);

real asinh(real x);
} // namespace seissol::dr::misc

#endif // SEISSOL_DR_MISC_H
