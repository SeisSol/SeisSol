#ifndef SEISSOL_SRC_SOLVER_MULTIPLESIMULATIONS_H_
#define SEISSOL_SRC_SOLVER_MULTIPLESIMULATIONS_H_

#include <cstddef>
#include <functional>
#include <init.h>
#include <yateto.h>

namespace seissol::multisim {

// duplicates the function argument `source` N times and calls `function` with it
template <std::size_t N, typename T, typename F, typename... Pack>
auto packed(F&& function, const T& source, Pack&&... copies) {
  if constexpr (sizeof...(Pack) < N) {
    return packed<N>(std::forward<F>(function), source, std::forward<Pack>(copies)..., source);
  } else {
    return std::invoke(std::forward<F>(function), std::forward<Pack>(copies)...);
  }
}

#ifdef MULTIPLE_SIMULATIONS
constexpr unsigned int NumSimulations = MULTIPLE_SIMULATIONS;
constexpr unsigned int BasisFunctionDimension = 1;
template <typename F, typename... Args>
auto multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), sim, std::forward<Args>(args)...);
}
template <typename T, typename F, typename... Args>
auto multisimObjectWrap(F&& func, T& obj, int sim, Args&&... args) {
  return std::invoke(std::forward<F>(func), obj, sim, std::forward<Args>(args)...);
}
template <unsigned Rank, typename RealT, typename IdxT>
auto simtensor(::yateto::DenseTensorView<Rank, RealT, IdxT>& tensor, int sim) {
  static_assert(Rank > 0, "Tensor rank needs to be non-scalar (rank > 0)");
  return packed<Rank - 1>([&](auto... args) { return tensor.subtensor(sim, args...); },
                          ::yateto::slice<>());
}
constexpr size_t MultisimStart = init::QAtPoint::Start[0];
constexpr size_t MultisimEnd = init::QAtPoint::Stop[0];
constexpr bool MultisimEnabled = true;

// last resort
#define SEISSOL_MULTISIM_WRAP(wrap, sim, ...) wrap(sim, __VA_ARGS__)
#else
constexpr unsigned int NumSimulations = 1;
constexpr unsigned int BasisFunctionDimension = 0;
template <typename F, typename... Args>
auto multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), std::forward<Args>(args)...);
}
template <typename T, typename F, typename... Args>
auto multisimObjectWrap(F&& func, T& obj, int sim, Args&&... args) {
  return std::invoke(std::forward<F>(func), obj, std::forward<Args>(args)...);
}
template <unsigned Rank, typename RealT, typename IdxT>
auto simtensor(::yateto::DenseTensorView<Rank, RealT, IdxT>& tensor, int sim) {
  return tensor;
}
constexpr size_t MultisimStart = 0;
constexpr size_t MultisimEnd = 1;
constexpr bool MultisimEnabled = false;

// last resort
#define SEISSOL_MULTISIM_WRAP(wrap, sim, ...) wrap(__VA_ARGS__)
#endif

} // namespace seissol::multisim

#endif // SEISSOL_SRC_SOLVER_MULTIPLESIMULATIONS_H_
