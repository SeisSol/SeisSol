// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Alignment.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include <cstdio>
#include <yateto.h>

#include <Solver/MultipleSimulations.h>

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

namespace {

__launch_bounds__(512) __global__ void kernel_local(const float** A,
                                                    const float** B,
                                                    unsigned Boffset,
                                                    const float* C1,
                                                    const float* C2,
                                                    const float* C3,
                                                    const float* C4,
                                                    float** D,
                                                    size_t numElements,
                                                    const unsigned* flags) {
  // meta data:
  // A = {rows: 64, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0), np.int64(64),
  // np.int64(9)]}; B = {rows: 9, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0),
  // np.int64(9), np.int64(9)]}; C = {rows: 56, cols: 56, addr: none, bbox: [np.int64(0),
  // np.int64(0), np.int64(56), np.int64(56)]}; D = {rows: 64, cols: 9, addr: pointer_based, bbox:
  // [np.int64(0), np.int64(0), np.int64(56), np.int64(9)]};

  // tmp0 = 1.0 * A x B
  // D = 1.0 * C x tmp0 + 1.0 * D

  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 64;
  constexpr int CountR = Count % 64;

  __shared__ __align__(8) float total_shrmem0[(576 + Count) * 8];

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;
  if (batchId < numElements) {
    const float* const __restrict__ glbA = A[batchId];
    const float* const __restrict__ glbB = B[batchId] + Boffset;
    const float* const __restrict__ glbC1 = C1;
    const float* const __restrict__ glbC2 = C2;
    const float* const __restrict__ glbC3 = C3;
    const float* const __restrict__ glbC4 = C4;
    float* const __restrict__ glbD = D[batchId];

    const bool has1 = (flags[batchId] & 1) != 0;
    const bool has2 = (flags[batchId] & 2) != 0;
    const bool has3 = (flags[batchId] & 4) != 0;
    const bool has4 = (flags[batchId] & 8) != 0;

    float reg0[Faces][Quantities]{};
    float reg1[Quantities]{};

    float* shrmem0 = &total_shrmem0[(576 + Count) * threadIdx.y];

    // writing to shr mem: from reg0 to _1
    float* __restrict__ _1 = &shrmem0[0];
#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      _1[tid_x + i * 64] = __builtin_nontemporal_load(&glbA[tid_x + i * 64]);
    }

    float* __restrict__ _0 = &shrmem0[576];

    // load all 4 9×9 matrices
    // loading glbB to _0: # no trans, extended
#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      _0[tid_x + i * 64] = glbB[tid_x + i * 64];
    }
    if (tid_x < CountR) {
      _0[tid_x + CountH * 64] = glbB[tid_x + CountH * 64];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      reg1[n] = __builtin_nontemporal_load(&glbD[tid_x + n * 64]);
    }

    // gemm: glbC x _1
    if (tid_x < 56) {

#pragma unroll 2
      for (int k = 0; k < 56; k += 8) {
        float values[Faces][8]{};
        if (has1) {
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
            values[0][kk] = glbC1[tid_x + (k + kk) * 56];
          }
        }
        if (has2) {
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
            values[1][kk] = glbC2[tid_x + (k + kk) * 56];
          }
        }
        if (has3) {
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
            values[2][kk] = glbC3[tid_x + (k + kk) * 56];
          }
        }
        if (has4) {
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
            values[3][kk] = glbC4[tid_x + (k + kk) * 56];
          }
        }

#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
          float local[8]{};
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
            local[kk] = _1[(k + kk) + n * 64];
          }
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
#pragma unroll
            for (int d = 0; d < Faces; ++d) {
              reg0[d][n] += values[d][kk] * local[kk];
            }
          }
        }
      }
    }

    // gemm: glbA x _0
#pragma unroll
    for (int d = 0; d < Faces; ++d) {
#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
#pragma unroll
        for (int k = 0; k < Quantities; ++k) {
          reg1[n] += reg0[d][k] * _0[k + n * Quantities + Quantities * Quantities * d];
        }
      }
    }

    // write results back to glb. memory
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg1[n], &glbD[tid_x + n * 64]);
    }
  }
}

} // namespace

namespace seissol::kernels::local::aux {
void launch_local(const float** A,
                  const float** B,
                  unsigned Boffset,
                  const float* C1,
                  const float* C2,
                  const float* C3,
                  const float* C4,
                  float** D,
                  size_t numElements,
                  const unsigned* flags,
                  void* streamPtr) {
  dim3 block(64, 8, 1);
  dim3 grid((numElements + 8 - 1) / 8, 1, 1);
  hipStream_t stream = (streamPtr != nullptr) ? static_cast<hipStream_t>(streamPtr) : 0;
  kernel_local<<<grid, block, 0, stream>>>(A, B, Boffset, C1, C2, C3, C4, D, numElements, flags);
}
} // namespace seissol::kernels::local::aux

namespace {

template <int Full, int Size1, int Size2>
__device__ __forceinline__ void kernel_ckstep(float reg1[9],
                                              float reg2[9],
                                              const float* C1,
                                              const float* C2,
                                              const float* C3,
                                              float** D,
                                              unsigned Doffset,
                                              float F,
                                              float* shrmem0) {
  // meta data:
  // A = {rows: 64, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0), np.int64(64),
  // np.int64(9)]}; B = {rows: 9, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0),
  // np.int64(9), np.int64(9)]}; C = {rows: 56, cols: 56, addr: none, bbox: [np.int64(0),
  // np.int64(0), np.int64(56), np.int64(56)]}; D = {rows: 64, cols: 9, addr: pointer_based, bbox:
  // [np.int64(0), np.int64(0), np.int64(56), np.int64(9)]};

  // tmp0 = 1.0 * A x B
  // D = 1.0 * C x tmp0 + 1.0 * D

  constexpr int Quantities = 9;
  constexpr int Faces = 3;

  constexpr int MyAlign = (seissol::Alignment / sizeof(float));
  constexpr int PaddedFull = ((Full + MyAlign - 1) / MyAlign) * MyAlign;
  constexpr int PaddedSize1 = ((Size1 + MyAlign - 1) / MyAlign) * MyAlign;
  constexpr int PaddedSize2 = ((Size2 + MyAlign - 1) / MyAlign) * MyAlign;

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;

  const float* const __restrict__ glbC1 = C1;
  const float* const __restrict__ glbC2 = C2;
  const float* const __restrict__ glbC3 = C3;
  float* const __restrict__ glbD = D[batchId] + Doffset;

  float reg0[Faces][Quantities]{};

  float* __restrict__ _0 = &shrmem0[PaddedSize1 * Quantities];
  float* __restrict__ _1 = &shrmem0[0];

  // gemm: glbC x _1
#pragma unroll 2
  for (int k = 0; k < (Size1 / 8) * 8; k += 8) {
    float values[Faces][8]{};
#pragma unroll
    for (int kk = 0; kk < 8; ++kk) {
      values[0][kk] = glbC1[tid_x + (k + kk) * PaddedFull];
    }
#pragma unroll
    for (int kk = 0; kk < 8; ++kk) {
      values[1][kk] = glbC2[tid_x + (k + kk) * PaddedFull];
    }
#pragma unroll
    for (int kk = 0; kk < 8; ++kk) {
      values[2][kk] = glbC3[tid_x + (k + kk) * PaddedFull];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      float local[8]{};
#pragma unroll
      for (int kk = 0; kk < 8; ++kk) {
        local[kk] = _1[(k + kk) + n * PaddedSize1];
      }
#pragma unroll
      for (int kk = 0; kk < 8; ++kk) {
#pragma unroll
        for (int d = 0; d < Faces; ++d) {
          reg0[d][n] += values[d][kk] * local[kk];
        }
      }
    }
  }

  if constexpr(Size1 % 8 != 0) {
    const int k = (Size1 / 8) * 8;
    float values[Faces][Size1 % 8]{};
#pragma unroll
    for (int kk = 0; kk < Size1 % 8; ++kk) {
      values[0][kk] = glbC1[tid_x + (k + kk) * PaddedFull];
    }
#pragma unroll
    for (int kk = 0; kk < Size1 % 8; ++kk) {
      values[1][kk] = glbC2[tid_x + (k + kk) * PaddedFull];
    }
#pragma unroll
    for (int kk = 0; kk < Size1 % 8; ++kk) {
      values[2][kk] = glbC3[tid_x + (k + kk) * PaddedFull];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      float local[Size1 % 8]{};
#pragma unroll
      for (int kk = 0; kk < Size1 % 8; ++kk) {
        local[kk] = _1[(k + kk) + n * PaddedSize1];
      }
#pragma unroll
      for (int kk = 0; kk < Size1 % 8; ++kk) {
#pragma unroll
        for (int d = 0; d < Faces; ++d) {
          reg0[d][n] += values[d][kk] * local[kk];
        }
      }
    }
  }

  // gemm: glbA x _0
#pragma unroll
  for (int d = 0; d < Faces; ++d) {
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
#pragma unroll
      for (int k = 0; k < Quantities; ++k) {
        reg1[n] += reg0[d][k] * _0[k + n * Quantities + Quantities * Quantities * d];
      }
    }
  }

#pragma unroll
  for (int n = 0; n < Quantities; ++n) {
    reg2[n] += F * reg1[n];
  }

  // write results back to glb. memory
  if (tid_x < PaddedSize2) {
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg1[n], &glbD[tid_x + n * PaddedSize2]);
    }
  }
}

template <int Full, int Size1>
__launch_bounds__(512) __global__ void kernel_cke(const float** A,
                                                  unsigned Aoffset,
                                                  const float** B,
                                                  unsigned Boffset,
                                                  const float* C1,
                                                  const float* C2,
                                                  const float* C3,
                                                  float** D,
                                                  unsigned Doffset,
                                                  float** E,
                                                  float F1,
                                                  float F2,
                                                  float F3,
                                                  size_t numElements) {
  // meta data:
  // A = {rows: 64, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0), np.int64(64),
  // np.int64(9)]}; B = {rows: 9, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0),
  // np.int64(9), np.int64(9)]}; C = {rows: 56, cols: 56, addr: none, bbox: [np.int64(0),
  // np.int64(0), np.int64(56), np.int64(56)]}; D = {rows: 64, cols: 9, addr: pointer_based, bbox:
  // [np.int64(0), np.int64(0), np.int64(56), np.int64(9)]};

  // tmp0 = 1.0 * A x B
  // D = 1.0 * C x tmp0 + 1.0 * D

  constexpr int Quantities = 9;
  constexpr int Faces = 3;

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 64;
  constexpr int CountR = Count % 64;

  constexpr int MyAlign = (seissol::Alignment / sizeof(float));
  constexpr int PaddedFull = ((Full + MyAlign - 1) / MyAlign) * MyAlign;
  constexpr int PaddedSize1 = ((Size1 + MyAlign - 1) / MyAlign) * MyAlign;

  __shared__ __align__(8) float total_shrmem0[(PaddedSize1 * Quantities + Count) * 8];

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;
  if (batchId < numElements) {
    const float* const __restrict__ glbA = A[batchId] + Aoffset;
    const float* const __restrict__ glbB = B[batchId] + Boffset;
    const float* const __restrict__ glbC1 = C1;
    const float* const __restrict__ glbC2 = C2;
    const float* const __restrict__ glbC3 = C3;
    float* const __restrict__ glbD = D[batchId] + Doffset;
    float* const __restrict__ glbE = E[batchId];

    float reg0[Faces][Quantities]{};
    float reg1[Quantities]{};
    float reg2[Quantities]{};

    float* shrmem0 = &total_shrmem0[(PaddedSize1 * Quantities + Count) * threadIdx.y];

    // writing to shr mem: from reg0 to _1
    float* __restrict__ _1 = &shrmem0[0];
    if (tid_x < PaddedSize1) {
#pragma unroll
      for (int i = 0; i < Quantities; ++i) {
        _1[tid_x + i * PaddedSize1] = __builtin_nontemporal_load(&glbA[tid_x + i * PaddedSize1]);
      }
    }

    float* __restrict__ _0 = &shrmem0[PaddedSize1 * Quantities];

    // load all 3 9×9 matrices
    // loading glbB to _0: # no trans, extended
#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      _0[tid_x + i * 64] = glbB[tid_x + i * 64];
    }
    if (tid_x < CountR) {
      _0[tid_x + CountH * 64] = glbB[tid_x + CountH * 64];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      reg2[n] = __builtin_nontemporal_load(&glbE[tid_x + n * PaddedFull]);
    }

    kernel_ckstep<Full, Size1, 10>(reg1, reg2, C1, C2, C3, D, Doffset, F1, shrmem0);
    if (tid_x < PaddedSize1) {
#pragma unroll
      for (int i = 0; i < Quantities; ++i) {
        _1[tid_x + i * PaddedSize1] = reg1[i];
      }
    }
    kernel_ckstep<Full, 10, 4>(reg1, reg2, C1, C2, C3, D, Doffset + 32 * 9, F2, shrmem0);
    if (tid_x < PaddedSize1) {
#pragma unroll
      for (int i = 0; i < Quantities; ++i) {
        _1[tid_x + i * PaddedSize1] = reg1[i];
      }
    }
    kernel_ckstep<Full, 4, 1>(reg1, reg2, C1, C2, C3, D, Doffset + 64 * 9, F3, shrmem0);

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg2[n], &glbE[tid_x + n * PaddedFull]);
    }
  }
}

template <int Full, int Size1, int Size2>
__launch_bounds__(512) __global__ void kernel_ck(const float** A,
                                                 unsigned Aoffset,
                                                 const float** B,
                                                 unsigned Boffset,
                                                 const float* C1,
                                                 const float* C2,
                                                 const float* C3,
                                                 float** D,
                                                 unsigned Doffset,
                                                 float** E,
                                                 float F,
                                                 size_t numElements) {
  // meta data:
  // A = {rows: 64, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0), np.int64(64),
  // np.int64(9)]}; B = {rows: 9, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0),
  // np.int64(9), np.int64(9)]}; C = {rows: 56, cols: 56, addr: none, bbox: [np.int64(0),
  // np.int64(0), np.int64(56), np.int64(56)]}; D = {rows: 64, cols: 9, addr: pointer_based, bbox:
  // [np.int64(0), np.int64(0), np.int64(56), np.int64(9)]};

  // tmp0 = 1.0 * A x B
  // D = 1.0 * C x tmp0 + 1.0 * D

  constexpr int Quantities = 9;
  constexpr int Faces = 3;

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 64;
  constexpr int CountR = Count % 64;

  constexpr int MyAlign = (seissol::Alignment / sizeof(float));
  constexpr int PaddedFull = ((Full + MyAlign - 1) / MyAlign) * MyAlign;
  constexpr int PaddedSize1 = ((Size1 + MyAlign - 1) / MyAlign) * MyAlign;
  constexpr int PaddedSize2 = ((Size2 + MyAlign - 1) / MyAlign) * MyAlign;

  __shared__ __align__(8) float total_shrmem0[(PaddedSize1 * Quantities + Count) * 8];

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;
  if (batchId < numElements) {
    const float* const __restrict__ glbA = A[batchId] + Aoffset;
    const float* const __restrict__ glbB = B[batchId] + Boffset;
    const float* const __restrict__ glbC1 = C1;
    const float* const __restrict__ glbC2 = C2;
    const float* const __restrict__ glbC3 = C3;
    float* const __restrict__ glbD = D[batchId] + Doffset;
    float* const __restrict__ glbE = E[batchId];

    float reg0[Faces][Quantities]{};
    float reg1[Quantities]{};
    float reg2[Quantities]{};

    float* shrmem0 = &total_shrmem0[(PaddedSize1 * Quantities + Count) * threadIdx.y];

    // writing to shr mem: from reg0 to _1
    float* __restrict__ _1 = &shrmem0[0];
    if (tid_x < PaddedSize1) {
#pragma unroll
      for (int i = 0; i < Quantities; ++i) {
        _1[tid_x + i * PaddedSize1] = __builtin_nontemporal_load(&glbA[tid_x + i * PaddedSize1]);
      }
    }

    float* __restrict__ _0 = &shrmem0[PaddedSize1 * Quantities];

    // load all 3 9×9 matrices
    // loading glbB to _0: # no trans, extended
#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      _0[tid_x + i * 64] = glbB[tid_x + i * 64];
    }
    if (tid_x < CountR) {
      _0[tid_x + CountH * 64] = glbB[tid_x + CountH * 64];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      reg2[n] = __builtin_nontemporal_load(&glbE[tid_x + n * PaddedFull]);
    }

    // gemm: glbC x _1
#pragma unroll 2
    for (int k = 0; k < (Size1 / 8) * 8; k += 8) {
      float values[Faces][8]{};
#pragma unroll
      for (int kk = 0; kk < 8; ++kk) {
        values[0][kk] = glbC1[tid_x + (k + kk) * PaddedFull];
      }
#pragma unroll
      for (int kk = 0; kk < 8; ++kk) {
        values[1][kk] = glbC2[tid_x + (k + kk) * PaddedFull];
      }
#pragma unroll
      for (int kk = 0; kk < 8; ++kk) {
        values[2][kk] = glbC3[tid_x + (k + kk) * PaddedFull];
      }

#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
        float local[8]{};
#pragma unroll
        for (int kk = 0; kk < 8; ++kk) {
          local[kk] = _1[(k + kk) + n * PaddedSize1];
        }
#pragma unroll
        for (int kk = 0; kk < 8; ++kk) {
#pragma unroll
          for (int d = 0; d < Faces; ++d) {
            reg0[d][n] += values[d][kk] * local[kk];
          }
        }
      }
    }

    if constexpr(Size1 % 8 != 0)
    {
      const int k = (Size1 / 8) * 8;
      float values[Faces][Size1 % 8]{};
#pragma unroll
      for (int kk = 0; kk < Size1 % 8; ++kk) {
        values[0][kk] = glbC1[tid_x + (k + kk) * PaddedFull];
      }
#pragma unroll
      for (int kk = 0; kk < Size1 % 8; ++kk) {
        values[1][kk] = glbC2[tid_x + (k + kk) * PaddedFull];
      }
#pragma unroll
      for (int kk = 0; kk < Size1 % 8; ++kk) {
        values[2][kk] = glbC3[tid_x + (k + kk) * PaddedFull];
      }

#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
        float local[Size1 % 8]{};
#pragma unroll
        for (int kk = 0; kk < Size1 % 8; ++kk) {
          local[kk] = _1[(k + kk) + n * PaddedSize1];
        }
#pragma unroll
        for (int kk = 0; kk < Size1 % 8; ++kk) {
#pragma unroll
          for (int d = 0; d < Faces; ++d) {
            reg0[d][n] += values[d][kk] * local[kk];
          }
        }
      }
    }

    // gemm: glbA x _0
#pragma unroll
    for (int d = 0; d < Faces; ++d) {
#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
#pragma unroll
        for (int k = 0; k < Quantities; ++k) {
          reg1[n] += reg0[d][k] * _0[k + n * Quantities + Quantities * Quantities * d];
        }
      }
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      reg2[n] += F * reg1[n];
    }

    // write results back to glb. memory
    if (tid_x < PaddedSize2) {
#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
        __builtin_nontemporal_store(reg1[n], &glbD[tid_x + n * PaddedSize2]);
      }
    }
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg2[n], &glbE[tid_x + n * PaddedFull]);
    }
  }
}

template <int Full>
__launch_bounds__(512) __global__
    void kernel_ick(const float** A, unsigned Aoffset, float** E, float F, size_t numElements) {
  // meta data:
  // A = {rows: 64, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0), np.int64(64),
  // np.int64(9)]}; B = {rows: 9, cols: 9, addr: pointer_based, bbox: [np.int64(0), np.int64(0),
  // np.int64(9), np.int64(9)]}; C = {rows: 56, cols: 56, addr: none, bbox: [np.int64(0),
  // np.int64(0), np.int64(56), np.int64(56)]}; D = {rows: 64, cols: 9, addr: pointer_based, bbox:
  // [np.int64(0), np.int64(0), np.int64(56), np.int64(9)]};

  // tmp0 = 1.0 * A x B
  // D = 1.0 * C x tmp0 + 1.0 * D

  constexpr int Quantities = 9;

  constexpr int MyAlign = (seissol::Alignment / sizeof(float));
  constexpr int PaddedFull = ((Full + MyAlign - 1) / MyAlign) * MyAlign;

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;
  if (batchId < numElements) {
    const float* const __restrict__ glbA = A[batchId] + Aoffset;
    float* const __restrict__ glbE = E[batchId];

    float reg1[Quantities]{};
    float reg2[Quantities]{};
#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      reg1[i] = __builtin_nontemporal_load(&glbA[tid_x + i * PaddedFull]);
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      reg2[n] = F * reg1[n];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg2[n], &glbE[tid_x + n * PaddedFull]);
    }
  }
}

} // namespace

namespace seissol::kernels::time::aux {

template <int SourceOrder>
void launch_cki(const float** A,
                unsigned Aoffset,
                const float** B,
                unsigned Boffset,
                const float* C1,
                const float* C2,
                const float* C3,
                float** D,
                unsigned Doffset,
                float** E,
                float F,
                size_t numElements,
                void* streamPtr) {
  dim3 block(64, 8, 1);
  dim3 grid((numElements + 8 - 1) / 8, 1, 1);
  hipStream_t stream = (streamPtr != nullptr) ? static_cast<hipStream_t>(streamPtr) : 0;
  constexpr std::size_t Full = seissol::kernels::getNumberOfBasisFunctions(ConvergenceOrder);
  constexpr std::size_t Size1 = seissol::kernels::getNumberOfBasisFunctions(SourceOrder);
  constexpr std::size_t Size2 = seissol::kernels::getNumberOfBasisFunctions(SourceOrder - 1);
  kernel_ck<Full, Size1, Size2><<<grid, block, 0, stream>>>(
      A, Aoffset, B, Boffset, C1, C2, C3, D, Doffset, E, F, numElements);
}

void launch_ck(const float** A,
               unsigned Aoffset,
               const float** B,
               unsigned Boffset,
               const float* C1,
               const float* C2,
               const float* C3,
               float** D,
               unsigned Doffset,
               float** E,
               float F,
               size_t numElements,
               void* streamPtr,
               int sourceOrder) {
  if (sourceOrder == 6) {
    launch_cki<6>(A, Aoffset, B, Boffset, C1, C2, C3, D, Doffset, E, F, numElements, streamPtr);
  }
  if (sourceOrder == 5) {
    launch_cki<5>(A, Aoffset, B, Boffset, C1, C2, C3, D, Doffset, E, F, numElements, streamPtr);
  }
  if (sourceOrder == 4) {
    launch_cki<4>(A, Aoffset, B, Boffset, C1, C2, C3, D, Doffset, E, F, numElements, streamPtr);
  }
  if (sourceOrder == 3) {
    launch_cki<3>(A, Aoffset, B, Boffset, C1, C2, C3, D, Doffset, E, F, numElements, streamPtr);
  }
  if (sourceOrder == 2) {
    launch_cki<2>(A, Aoffset, B, Boffset, C1, C2, C3, D, Doffset, E, F, numElements, streamPtr);
  }
}

void launch_ick(
    const float** A, unsigned Aoffset, float** E, float F, size_t numElements, void* streamPtr) {
  dim3 block(64, 8, 1);
  dim3 grid((numElements + 8 - 1) / 8, 1, 1);
  hipStream_t stream = (streamPtr != nullptr) ? static_cast<hipStream_t>(streamPtr) : 0;
  constexpr std::size_t Full = seissol::kernels::getNumberOfBasisFunctions(ConvergenceOrder);
  kernel_ick<Full><<<grid, block, 0, stream>>>(A, Aoffset, E, F, numElements);
}

} // namespace seissol::kernels::time::aux

#endif

namespace {
using namespace seissol::multisim;

template <typename Tensor>
constexpr size_t leadDim() {
  if constexpr (MultisimEnabled) {
    return Tensor::Stop[1] - Tensor::Start[1];
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Tensor>
constexpr size_t linearDim() {
  if constexpr (MultisimEnabled) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

constexpr auto getblock(int size) {
  if constexpr (MultisimEnabled) {
    return dim3(NumSimulations, size);
  } else {
    return dim3(size);
  }
}

__forceinline__ __device__ auto linearidx() {
  if constexpr (MultisimEnabled) {
    return threadIdx.y * NumSimulations + threadIdx.x;
  } else {
    return threadIdx.x;
  }
}

__forceinline__ __device__ auto linearsize() {
  if constexpr (MultisimEnabled) {
    return blockDim.y * blockDim.x;
  } else {
    return blockDim.x;
  }
}

__forceinline__ __device__ auto simidx() {
  if constexpr (MultisimEnabled) {
    return threadIdx.x;
  } else {
    return 0;
  }
}

__forceinline__ __device__ auto validx() {
  if constexpr (MultisimEnabled) {
    return threadIdx.y;
  } else {
    return threadIdx.x;
  }
}
} // namespace

namespace seissol::kernels::local_flux::aux::details {

__global__ void kernelFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                                         real** displacementDataPtrs,
                                         double* rhos,
                                         double g,
                                         size_t numElements) {

  const int tid = linearidx();
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    const double rho = rhos[elementId];
    real* elementBoundaryDofs = dofsFaceBoundaryNodalPtrs[elementId];
    real* elementDisplacement = displacementDataPtrs[elementId];

    constexpr auto numNodes = linearDim<seissol::nodal::init::nodes2D>();
    if (tid < numNodes) {
      constexpr auto ldINodal = linearDim<seissol::init::INodal>();

      const auto pressureAtBnd = static_cast<real>(-1.0) * rho * g * elementDisplacement[tid];

#pragma unroll
      for (int component{0}; component < 3; ++component) {
        elementBoundaryDofs[tid + component * ldINodal] =
            2.0 * pressureAtBnd - elementBoundaryDofs[tid + component * ldINodal];
      }
    }
  }
}

void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream) {
  dim3 block = getblock(leadDim<seissol::nodal::init::nodes2D>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelFreeSurfaceGravity<<<grid, block, 0, stream>>>(
      dofsFaceBoundaryNodalPtrs, displacementDataPtrs, rhos, g, numElements);
}

__global__ void kernelEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                                   real** easiBoundaryMapPtrs,
                                   real** easiBoundaryConstantPtrs,
                                   size_t numElements) {

  const int tid = validx();
  const int elementId = blockIdx.x;

  constexpr auto ldINodalDim = linearDim<seissol::init::INodal>();
  constexpr auto iNodalDim0 = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension + 0];
  constexpr auto iNodalDim1 = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension + 1];
  __shared__ __align__(8) real resultTerm[iNodalDim1][iNodalDim0][multisim::NumSimulations];

  constexpr auto ldConstantDim = linearDim<seissol::init::easiBoundaryConstant>();
  constexpr auto constantDim0 =
      seissol::tensor::easiBoundaryConstant::Shape[multisim::BasisFunctionDimension + 0];
  constexpr auto constantDim1 =
      seissol::tensor::easiBoundaryConstant::Shape[multisim::BasisFunctionDimension + 1];
  __shared__ __align__(8) real rightTerm[constantDim0][constantDim1][multisim::NumSimulations];

  constexpr auto ldMapDim = leadDim<seissol::init::easiBoundaryMap>();
  constexpr auto mapDim0 = seissol::tensor::easiBoundaryMap::Shape[0];
  constexpr auto mapDim1 = seissol::tensor::easiBoundaryMap::Shape[1];
  constexpr auto mapDim2 = seissol::tensor::easiBoundaryMap::Shape[2];
  __shared__ __align__(8) real leftTerm[mapDim0][mapDim2];

  static_assert(iNodalDim1 == constantDim0, "supposed to be equal");
  static_assert(iNodalDim1 == mapDim0, "supposed to be equal");

  if (elementId < numElements) {
    real* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
    real* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
    auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

    for (int i = linearidx(); i < (ldConstantDim * constantDim1); i += linearsize()) {
      const auto sim = i % multisim::NumSimulations;
      const auto subsim = i / multisim::NumSimulations;
      const auto quantity = subsim % constantDim0;
      const auto quadpoint = subsim / constantDim0;
      rightTerm[quantity][quadpoint][sim] = easiBoundaryConstant[i];
    }
    __syncthreads();

    for (int i = 0; i < iNodalDim1; ++i) {
      if (tid < iNodalDim0) {
        resultTerm[i][tid][simidx()] = 0.0;
      }
    }
    __syncthreads();

    for (int quantity = 0; quantity < mapDim1; ++quantity) {
      for (int quadpoint = 0; quadpoint < mapDim2; ++quadpoint) {
        if (tid < mapDim0) {
          leftTerm[tid][quadpoint] =
              easiBoundaryMap[tid + ldMapDim * (quantity + quadpoint * mapDim1)];
        }
      }
      __syncthreads();

      if (tid < mapDim2) {
        const real col = dofsFaceBoundaryNodal[linearidx() + quantity * ldINodalDim];
        for (int quantity2 = 0; quantity2 < mapDim0; ++quantity2) {
          resultTerm[quantity2][tid][simidx()] += leftTerm[quantity2][tid] * col;
        }
      }
      __syncthreads();
    }

    if (tid < iNodalDim0) {
      for (int quantity2 = 0; quantity2 < iNodalDim1; ++quantity2) {
        dofsFaceBoundaryNodal[linearidx() + quantity2 * ldINodalDim] =
            resultTerm[quantity2][tid][simidx()] + rightTerm[quantity2][tid][simidx()];
      }
    }
  }
}

void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  dim3 block = getblock(leadDim<seissol::init::INodal>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelEasiBoundary<<<grid, block, 0, stream>>>(
      dofsFaceBoundaryNodalPtrs, easiBoundaryMapPtrs, easiBoundaryConstantPtrs, numElements);
}
} // namespace seissol::kernels::local_flux::aux::details

namespace seissol::kernels::time::aux {
__global__ void kernelextractRotationMatrices(real** displacementToFaceNormalPtrs,
                                              real** displacementToGlobalDataPtrs,
                                              real** TPtrs,
                                              real** TinvPtrs,
                                              size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    real* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
    real* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
    auto* T = TPtrs[elementId];
    auto* Tinv = TinvPtrs[elementId];

    constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv>();
    constexpr auto ldT = yateto::leadDim<seissol::init::T>();
    constexpr auto ldDisplacement = yateto::leadDim<seissol::init::displacementRotationMatrix>();

    const int i = threadIdx.x;
    const int j = threadIdx.y;

    displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
    displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
  }
}

void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream) {
  dim3 block(3, 3, 1);
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelextractRotationMatrices<<<grid, block, 0, stream>>>(
      displacementToFaceNormalPtrs, displacementToGlobalDataPtrs, TPtrs, TinvPtrs, numElements);
}

__global__ void
    kernelInitializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                         real** integratedDisplacementNodalPtrs,
                                                         real** rotatedFaceDisplacementPtrs,
                                                         double deltaTInt,
                                                         size_t numElements) {

  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    auto* prevCoefficients = prevCoefficientsPtrs[elementId];
    auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
    const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

    assert(linearDim<seissol::nodal::init::nodes2D>() <=
           linearDim<seissol::init::rotatedFaceDisplacement>());

    const int tid = linearidx();
    constexpr auto num2dNodes = linearDim<seissol::nodal::init::nodes2D>();
    if (tid < num2dNodes) {
      prevCoefficients[tid] = rotatedFaceDisplacement[tid];
      integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
    }
  }
}

void initializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** rotatedFaceDisplacementPtrs,
                                                    double deltaTInt,
                                                    size_t numElements,
                                                    void* deviceStream) {

  dim3 block = getblock(leadDim<seissol::nodal::init::nodes2D>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelInitializeTaylorSeriesForGravitationalBoundary<<<grid, block, 0, stream>>>(
      prevCoefficientsPtrs,
      integratedDisplacementNodalPtrs,
      rotatedFaceDisplacementPtrs,
      deltaTInt,
      numElements);
}

__global__ void kernelComputeInvAcousticImpedance(double* invImpedances,
                                                  double* rhos,
                                                  double* lambdas,
                                                  size_t numElements) {

  size_t index = threadIdx.x + blockIdx.x * blockDim.x;
  if (index < numElements) {
    invImpedances[index] = 1.0 / std::sqrt(lambdas[index] * rhos[index]);
  }
}

void computeInvAcousticImpedance(
    double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream) {
  constexpr size_t blockSize{1024};
  dim3 block(blockSize, 1, 1);
  dim3 grid((numElements + blockSize - 1) / blockSize, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelComputeInvAcousticImpedance<<<grid, block, 0, stream>>>(
      invImpedances, rhos, lambdas, numElements);
}

__global__ void kernelUpdateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                                    real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** dofsFaceNodalPtrs,
                                                    double* invImpedances,
                                                    double* rhos,
                                                    double g,
                                                    double factorEvaluated,
                                                    double factorInt,
                                                    size_t numElements) {
  const int elementId = blockIdx.x;
  if (elementId < numElements) {
    constexpr int pIdx = 0;
    constexpr int uIdx = model::MaterialT::TractionQuantities;
    constexpr auto num2dNodes = linearDim<seissol::nodal::init::nodes2D>();

    const int tid = linearidx();
    if (tid < num2dNodes) {

      real* dofsFaceNodal = dofsFaceNodalPtrs[elementId];
      constexpr auto ldINodal = linearDim<seissol::init::INodal>();

      const auto uInside = dofsFaceNodal[tid + (uIdx + 0) * ldINodal];
      const auto vInside = dofsFaceNodal[tid + (uIdx + 1) * ldINodal];
      const auto wInside = dofsFaceNodal[tid + (uIdx + 2) * ldINodal];
      const auto pressureInside = dofsFaceNodal[tid + pIdx * ldINodal];

      real* prevCoefficients = prevCoefficientsPtrs[elementId];
      const auto rho = rhos[elementId];
      const auto invImpedance = invImpedances[elementId];

      const double curCoeff =
          uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
      prevCoefficients[tid] = curCoeff;

      constexpr auto ldFaceDisplacement = linearDim<seissol::init::faceDisplacement>();
      static_assert(num2dNodes <= ldFaceDisplacement, "");

      real* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];
      rotatedFaceDisplacement[tid + 0 * ldFaceDisplacement] += factorEvaluated * curCoeff;
      rotatedFaceDisplacement[tid + 1 * ldFaceDisplacement] += factorEvaluated * vInside;
      rotatedFaceDisplacement[tid + 2 * ldFaceDisplacement] += factorEvaluated * wInside;

      constexpr auto ldIntegratedFaceDisplacement =
          linearDim<seissol::init::averageNormalDisplacement>();
      static_assert(num2dNodes <= ldIntegratedFaceDisplacement, "");

      real* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
      integratedDisplacementNodal[tid] += factorInt * curCoeff;
    }
  }
}

void updateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                   real** prevCoefficientsPtrs,
                                   real** integratedDisplacementNodalPtrs,
                                   real** dofsFaceNodalPtrs,
                                   double* invImpedances,
                                   double* rhos,
                                   double g,
                                   double factorEvaluated,
                                   double factorInt,
                                   size_t numElements,
                                   void* deviceStream) {
  dim3 block = getblock(leadDim<seissol::nodal::init::nodes2D>());
  dim3 grid(numElements, 1, 1);
  auto stream = reinterpret_cast<StreamT>(deviceStream);
  kernelUpdateRotatedFaceDisplacement<<<grid, block, 0, stream>>>(rotatedFaceDisplacementPtrs,
                                                                  prevCoefficientsPtrs,
                                                                  integratedDisplacementNodalPtrs,
                                                                  dofsFaceNodalPtrs,
                                                                  invImpedances,
                                                                  rhos,
                                                                  g,
                                                                  factorEvaluated,
                                                                  factorInt,
                                                                  numElements);
}
} // namespace seissol::kernels::time::aux
