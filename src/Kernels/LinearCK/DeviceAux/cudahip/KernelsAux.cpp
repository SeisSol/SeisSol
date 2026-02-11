// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Solver/MultipleSimulations.h"

#include <cstdio>
#include <yateto.h>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

#ifdef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
namespace {
constexpr std::size_t Blocksize = 128;

template <typename SourceRealT>
constexpr std::size_t RestFunctions(std::size_t SourceOrder, std::size_t ThisOrder) {
  std::size_t total = 0;
  for (std::size_t j = ThisOrder; j < SourceOrder; ++j) {
    total +=
        seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  }
  return total;
}

template <std::size_t Quantities,
          std::size_t ThisOrder,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t Offset,
          std::size_t SharedOffset,
          typename... Coeffs>
static __device__ __forceinline__ void taylorSumInner(TargetRealT* const __restrict__ target,
                                                      const SourceRealT* const __restrict__ source,
                                                      SourceRealT* const __restrict__ shmem,
                                                      TargetRealT reg[Quantities],
                                                      TargetRealT coeff,
                                                      Coeffs... coeffs) {
  constexpr std::size_t MemorySize =
      seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) * Quantities;
  constexpr std::size_t SourceStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);
  constexpr std::size_t RestMemSize =
      SharedOffset > 0 ? 0 : RestFunctions<SourceRealT>(SourceOrder, ThisOrder) * Quantities;
  constexpr std::size_t SourceMemSize = SourceStride * Quantities;
  constexpr bool UseShared = MemorySize >= RestMemSize;
  constexpr std::size_t LoadSize = UseShared ? RestMemSize : SourceMemSize;

  static_assert(seissol::tensor::dQ::size(ThisOrder) == SourceMemSize,
                "Tensor size mismatch in explicit kernel.");

  if constexpr (LoadSize > 0) {
    __syncthreads();
    constexpr std::size_t Rounds = LoadSize / Blocksize;
    constexpr std::size_t Rest = LoadSize % Blocksize;
    if constexpr (Rounds > 0) {
#pragma unroll
      for (std::size_t j = 0; j < Rounds; ++j) {
        shmem[j * Blocksize + threadIdx.x] = source[Offset + j * Blocksize + threadIdx.x];
      }
    }
    if constexpr (Rest > 0) {
      if (threadIdx.x < Rest) {
        shmem[Rounds * Blocksize + threadIdx.x] = source[Offset + Rounds * Blocksize + threadIdx.x];
      }
    }
    __syncthreads();
  }

  constexpr std::size_t BasisFunctionsSize = std::min(SourceStride, TargetStride);
  if (threadIdx.x < BasisFunctionsSize) {
#pragma unroll
    for (std::size_t j = 0; j < Quantities; ++j) {
      // FIXME: non-optimal warp utilization (as before... But now it's in registers)
      reg[j] +=
          coeff * static_cast<TargetRealT>(shmem[SharedOffset + SourceStride * j + threadIdx.x]);
    }
  }
  // TODO(David): are we sure about this? Or is ThisOrder + 1 < SourceOrder enough?
  if constexpr (ThisOrder + 1 < std::min(SourceOrder, TargetOrder) && sizeof...(Coeffs) > 0) {
    constexpr std::size_t SharedPosition = UseShared ? SharedOffset + SourceMemSize : 0;
    taylorSumInner<Quantities,
                   ThisOrder + 1,
                   SourceRealT,
                   TargetRealT,
                   SourceOrder,
                   TargetOrder,
                   Offset + LoadSize,
                   SharedPosition>(target, source, shmem, reg, coeffs...);
  }
}

template <std::size_t Quantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          typename... Coeffs>
static __global__ __launch_bounds__(Blocksize) void taylorSumKernel(TargetRealT** targetBatch,
                                                                    const SourceRealT** sourceBatch,
                                                                    Coeffs... coeffs) {
  int batchId = blockIdx.x;

  __shared__ SourceRealT
      shmem[seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) *
            Quantities];
  TargetRealT reg[Quantities] = {0};

  const SourceRealT* const __restrict__ source =
      const_cast<const SourceRealT*>(sourceBatch[batchId]);
  TargetRealT* const __restrict__ target = targetBatch[batchId];

  taylorSumInner<Quantities, 0, SourceRealT, TargetRealT, SourceOrder, TargetOrder, 0, 0>(
      target, source, shmem, reg, coeffs...);

  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  if (threadIdx.x < TargetStride) {
#pragma unroll
    for (std::size_t j = 0; j < Quantities; ++j) {
      target[TargetStride * j + threadIdx.x] = reg[j];
    }
  }
}

template <std::size_t SourceQuantities,
          std::size_t TargetQuantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t... Indices>
void static taylorSumInternal(std::size_t count,
                              TargetRealT** target,
                              const SourceRealT** source,
                              const TargetRealT* coeffs,
                              void* stream,
                              std::index_sequence<Indices...> indices) {
  constexpr std::size_t Quantities = std::min(SourceQuantities, TargetQuantities);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  dim3 threads(Blocksize);
  dim3 blocks(count);

  StreamT castedStream = reinterpret_cast<StreamT>(stream);

  taylorSumKernel<Quantities, SourceRealT, TargetRealT, SourceOrder, TargetOrder>
      <<<blocks, threads, 0, castedStream>>>(target, source, coeffs[Indices]...);
}
} // namespace

namespace seissol::kernels::time::aux {
void taylorSum(
    std::size_t count, real** target, const real** source, const real* coeffs, void* stream) {
  taylorSumInternal<seissol::model::MaterialT::NumQuantities,
                    seissol::model::MaterialT::NumQuantities,
                    real,
                    real,
                    ConvergenceOrder,
                    ConvergenceOrder>(
      count, target, source, coeffs, stream, std::make_index_sequence<ConvergenceOrder>());
}
} // namespace seissol::kernels::time::aux
#endif
