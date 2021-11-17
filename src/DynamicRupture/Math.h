#ifndef SEISSOL_MATH_H
#define SEISSOL_MATH_H

#define DR_AS_QUADRATURE 0
#define DR_AS_CELLAVERAGE 1

#include <stdexcept>
#include <cmath>

namespace seissol::dr::aux {
template <int DrMethod>
constexpr int numGaussPoints2d(const int Order) {
  throw std::runtime_error("unknown Dynamic Rupture method");
  return -1;
}

template <>
constexpr int numGaussPoints2d<DR_AS_QUADRATURE>(const int Order) {
  return std::pow(Order + 1, 2);
}

struct CellAveragePrecomputed {
  // Precomputed values according to:
  // numberOfPoints = int(4**math.ceil(math.log(order*(order+1)/2,4)))
  constexpr static int NumberOfPoints[10] = {1, 4, 16, 16, 16, 64, 64, 64, 64, 64};

  constexpr static int getNumPoints(int Order) { return NumberOfPoints[Order]; }
};

template <>
constexpr int numGaussPoints2d<DR_AS_CELLAVERAGE>(const int Order) {
  return CellAveragePrecomputed::getNumPoints(Order);
}

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
} // namespace seissol::dr::aux

#endif // SEISSOL_MATH_H
