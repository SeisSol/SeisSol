// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_MULTIPLESIMULATIONS_H_
#define SEISSOL_SRC_SOLVER_MULTIPLESIMULATIONS_H_

#include "Config.h"
#include "GeneratedCode/init.h"

#include <cstddef>
#include <functional>
#include <tuple>
#include <utility>
#include <yateto.h>

// disable the omp simd declarations for old Intel compilers

#ifdef __INTEL_COMPILER
#define SEISSOL_NO_OMPSIMD
#endif // __INTEL_COMPILER

#ifdef __INTEL_LLVM_COMPILER
#if __INTEL_LLVM_COMPILER < 20230000
#define SEISSOL_NO_OMPSIMD
#endif
#endif // __INTEL_LLVM_COMPILER

namespace seissol::multisim {

// duplicates the function argument `source` N times and calls `function` with it
template <std::size_t N, typename T, typename F, typename... Pack>
decltype(auto) packed(F&& function, const T& source, Pack&&... copies) {
  if constexpr (sizeof...(Pack) < N) {
    return packed<N>(std::forward<F>(function), source, std::forward<Pack>(copies)..., source);
  } else {
    return std::invoke(std::forward<F>(function), std::forward<Pack>(copies)...);
  }
}

template <std::size_t Idx, typename F, typename T1, typename T2>
decltype(auto) reverseCallInternal(F&& function, T1&& fwdTuple, T2&& bckTuple) {
  if constexpr (std::tuple_size_v<T1> == Idx) {
    return std::apply(std::forward<F>(function), std::forward<T2>(bckTuple));
  } else {
    auto newBck =
        std::tuple_cat(std::make_tuple(std::get<Idx>(fwdTuple)), std::forward<T2>(bckTuple));
    return reverseCallInternal<Idx + 1>(
        std::forward<F>(function), std::forward<T1>(fwdTuple), newBck);
  }
}

// reverse the parameter pack arguments (`values`) and call `function` with the reversed pack
template <typename F, typename... Pack>
decltype(auto) reverseCall(F&& function, Pack&&... values) {
  std::tuple<> emptytuple{};
  return reverseCallInternal<0>(
      std::forward<F>(function), std::forward_as_tuple(std::forward<Pack>(values)...), emptytuple);
}

/// The helpers for the fused simulations of the configuration `Cfg`. `NumSimulationsT` only selects
/// the case and keeps its default.
template <typename Cfg, unsigned int NumSimulationsT = Cfg::NumSimulations>
struct MultisimHelperWrapper {
  // the (non-?)default case: NumSimulations > 1
  constexpr static unsigned int NumSimulations = NumSimulationsT;
  constexpr static unsigned int BasisFunctionDimension = 1;

  // The simulation index is the leading dimension of the fused tensors, and the hand-written parts
  // of SeisSol step through it with NumSimulations as the stride. So the code generator must not
  // pad it; codegen/generate.py chooses the vector size accordingly.
  static_assert(init::Q<Cfg>::Stop[0] - init::Q<Cfg>::Start[0] == NumSimulationsT,
                "The simulation dimension of the fused tensors is padded. Choose a vector size "
                "that divides the fused simulations (in bytes).");

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
  template <typename F, typename... Args>
  static decltype(auto) multisimWrap(F&& function, size_t sim, Args&&... args) {
    return std::invoke(std::forward<F>(function), sim, std::forward<Args>(args)...);
  }

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
  template <typename T, typename F, typename... Args>
  static decltype(auto) multisimObjectWrap(F&& func, T& obj, int sim, Args&&... args) {
    return std::invoke(std::forward<F>(func), obj, sim, std::forward<Args>(args)...);
  }

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
  template <typename F, typename... Args>
  static decltype(auto) multisimTranspose(F&& function, Args&&... args) {
    return reverseCall(std::forward<F>(function), std::forward<Args>(args)...);
  }

  template <typename TensorViewT>
  static decltype(auto) simtensor(TensorViewT& tensor, int sim) {
    static_assert(TensorViewT::dim() > 0, "Tensor rank needs to be non-scalar (rank > 0)");
    return packed<TensorViewT::dim() - 1>(
        [&](auto... args) { return tensor.subtensor(sim, args...); }, ::yateto::slice<>());
  }

  constexpr static size_t MultisimStart = init::QAtPoint<Cfg>::Start[0];
  constexpr static size_t MultisimEnd = init::QAtPoint<Cfg>::Stop[0];
  constexpr static bool MultisimEnabled = true;
};

template <typename Cfg>
struct MultisimHelperWrapper<Cfg, 1> {
  constexpr static unsigned int NumSimulations = 1;
  constexpr static unsigned int BasisFunctionDimension = 0;

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
  template <typename F, typename... Args>
  static decltype(auto) multisimWrap(F&& function, size_t /*sim*/, Args&&... args) {
    return std::invoke(std::forward<F>(function), std::forward<Args>(args)...);
  }

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
  template <typename T, typename F, typename... Args>
  static decltype(auto) multisimObjectWrap(F&& func, T& obj, int /*sim*/, Args&&... args) {
    return std::invoke(std::forward<F>(func), obj, std::forward<Args>(args)...);
  }

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
  template <typename F, typename... Args>
  static decltype(auto) multisimTranspose(F&& function, Args&&... args) {
    return std::invoke(std::forward<F>(function), std::forward<Args>(args)...);
  }

  template <typename TensorViewT>
  static decltype(auto) simtensor(TensorViewT& tensor, int /*sim*/) {
    return tensor;
  }
  constexpr static size_t MultisimStart = 0;
  constexpr static size_t MultisimEnd = 1;
  constexpr static bool MultisimEnabled = false;
};

/// The dimension of the basis functions in the tensors of the configuration `Cfg`.
template <typename Cfg>
constexpr unsigned int BasisDim = MultisimHelperWrapper<Cfg>::BasisFunctionDimension;

// The functions below take the configuration whose fused simulations they work on.

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
template <typename Cfg, typename F, typename... Args>
decltype(auto) multisimWrap(F&& function, size_t sim, Args&&... args) {
  return MultisimHelperWrapper<Cfg>::multisimWrap(
      std::forward<F>(function), sim, std::forward<Args>(args)...);
}

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
template <typename Cfg, typename T, typename F, typename... Args>
decltype(auto) multisimObjectWrap(F&& func, T& obj, int sim, Args&&... args) {
  return MultisimHelperWrapper<Cfg>::multisimObjectWrap(
      std::forward<F>(func), obj, sim, std::forward<Args>(args)...);
}

#ifndef SEISSOL_NO_OMPSIMD
#pragma omp declare simd
#endif
template <typename Cfg, typename F, typename... Args>
decltype(auto) multisimTranspose(F&& function, Args&&... args) {
  return MultisimHelperWrapper<Cfg>::multisimTranspose(std::forward<F>(function),
                                                       std::forward<Args>(args)...);
}

template <typename Cfg, typename TensorViewT>
decltype(auto) simtensor(TensorViewT& tensor, int sim) {
  return MultisimHelperWrapper<Cfg>::simtensor(tensor, sim);
}
template <typename Cfg, typename Tensor>
constexpr size_t leadDim() {
  if constexpr (MultisimHelperWrapper<Cfg>::MultisimEnabled) {
    return Tensor::Stop[1] - Tensor::Start[1];
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Cfg, typename Tensor>
constexpr size_t linearDim() {
  if constexpr (MultisimHelperWrapper<Cfg>::MultisimEnabled) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

} // namespace seissol::multisim

#endif // SEISSOL_SRC_SOLVER_MULTIPLESIMULATIONS_H_
