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
// TODO for FP64: use DPP64 maybe?
// i.e. we can broadcast per 16 threads
// however, no idea how to stuff constant matrices...
// (we'd need 128 KiB LDS for that at least)

template <typename T>
__device__ __forceinline__ auto readlane(T value, int lane) -> T {
  static_assert(sizeof(T) == sizeof(int), "NYI");
  int vc = *reinterpret_cast<int*>(&value);
  int vcr = __builtin_amdgcn_readlane(vc, lane);
  return *reinterpret_cast<T*>(&vcr);
}

template <int dpp1, int dpp2, int dpp3, bool dpp4, typename T>
__device__ __forceinline__ auto dpp(T value) -> T {
  static_assert(sizeof(T) == sizeof(int), "NYI");
  int vc = *reinterpret_cast<int*>(&value);
  int vcr = __builtin_amdgcn_mov_dpp(vc, dpp1, dpp2, dpp3, dpp4);
  return *reinterpret_cast<T*>(&vcr);
}

template <int lane, typename T>
__device__ __forceinline__ auto dpp4(T value) -> T {
  return dpp<(lane << 6) | (lane << 4) | (lane << 2) | lane, 0b1111, 0b1111, true>(value);
}

using af4 = __attribute__((__vector_size__(4 * sizeof(float)))) float;

#define ISTRINGIFY(x) #x
#define STR(x) ISTRINGIFY(x)
#define CM4STR(p1, p2, p3, p4, c, a, b)                                                            \
  "v_cndmask_b32_dpp " c ", " a ", " b ", vcc quad_perm:[" STR(p1) "," STR(p2) "," STR(            \
      p3) "," STR(p4) "] row_mask:0xf bank_mask:0xf bound_ctrl:1"
#define CMRSTR(cnt, c, a, b)                                                                       \
  "v_cndmask_b32_dpp " c ", " a ", " b                                                             \
  ", vcc row_ror:" STR(cnt) " row_mask:0xf bank_mask:0xf bound_ctrl:1"

template <typename T>
__device__ __forceinline__ auto
    transpose4x4bcst(T v1, T v2, T v3, T v4, int i1, int i2, int i3, int i4) -> T {
  const auto r1 = readlane(v1, i1);
  const auto r2 = readlane(v2, i2);
  const auto r3 = readlane(v3, i3);
  const auto r4 = readlane(v4, i4);

  const uint64_t mask1a = 0x1111111111111111ULL;
  const uint64_t mask1b = 0x4444444444444444ULL;
  const uint64_t mask2 = 0x3333333333333333ULL;

  T w1;
  T w2;
  T w;

  // NOTE: we upcast SRC0 AND SRC1 to VGPR, because of constant bus restrictions
  // (on e.g. gfx1103 wave32, one of SRC0/SRC1 could stay SGPR)

  __asm("v_cndmask_b32_e64 %[w], %[r2], %[r1], %[mask]"
        : [w] "=v"(w1)
        : [r2] "v"(r2), [r1] "v"(r1), [mask] "s"(mask1a)
        :);
  __asm("v_cndmask_b32_e64 %[w], %[r2], %[r1], %[mask]"
        : [w] "=v"(w2)
        : [r2] "v"(r4), [r1] "v"(r3), [mask] "s"(mask1b)
        :);
  __asm("v_cndmask_b32_e64 %[w], %[r2], %[r1], %[mask]"
        : [w] "=v"(w)
        : [r2] "v"(w2), [r1] "v"(w1), [mask] "s"(mask2)
        :);

  return w;
}

template <typename T>
__device__ __forceinline__ auto
    transpose16x4b4x1(T& w1, T& w2, T& w3, T& w4, T v1, T v2, T v3, T v4) {

  // TODO: rewrite with DPP movs (possible here via row modifiers)

  const uint64_t mask1a = 0x0f0f0f0f0f0f0f0fULL;
  const uint64_t mask1b = 0xf0f0f0f0f0f0f0f0ULL;
  const uint64_t mask2a = 0x00ff00ff00ff00ffULL;
  const uint64_t mask2b = 0xff00ff00ff00ff00ULL;

  T u1, u2, u3, u4;

  // 11 12 13 14
  // 21 22 23 24
  // 31 32 33 34
  // 41 42 43 44

  // 11 21 13 23 (DPP for row 2)
  // 12 22 14 24 (DPP for row 1)
  // 31 41 33 43 (DPP for row 4)
  // 32 42 34 44 (DPP for row 3)

  // 11 21 31 41 (DPP for row 3)
  // 12 22 32 42 (DPP for row 4)
  // 13 23 33 43 (DPP for row 1)
  // 14 24 34 44 (DPP for row 2)

  // dpp<0x120 + offset, 0b0101, 0b1111, false>();

  // clang-format off

  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CMRSTR(12, "%[u1]", "%[v2]", "%[v1]") "\n\t"
  CMRSTR(12, "%[u3]", "%[v4]", "%[v3]") "\n\t"
  : [u1] "=v" (u1), [u3] "=v" (u3)
  : [mask] "s" (mask1a), [v1] "v" (v1), [v2] "v" (v2), [v3] "v" (v3), [v4] "v" (v4)
  : "vcc");
  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CMRSTR(4, "%[u2]", "%[v1]", "%[v2]") "\n\t"
  CMRSTR(4, "%[u4]", "%[v3]", "%[v4]") "\n\t"
  : [u2] "=v" (u2), [u4] "=v" (u4)
  : [mask] "s" (mask1b), [v1] "v" (v1), [v2] "v" (v2), [v3] "v" (v3), [v4] "v" (v4)
  : "vcc");
  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CMRSTR(8, "%[w1]", "%[u3]", "%[u1]") "\n\t"
  CMRSTR(8, "%[w2]", "%[u4]", "%[u2]") "\n\t"
  : [w1] "=v" (w1), [w2] "=v" (w2)
  : [mask] "s" (mask2a), [u1] "v" (u1), [u2] "v" (u2), [u3] "v" (u3), [u4] "v" (u4)
  : "vcc");
  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CMRSTR(8, "%[w3]", "%[u1]", "%[u3]") "\n\t"
  CMRSTR(8, "%[w4]", "%[u2]", "%[u4]")
  : [w3] "=v" (w3), [w4] "=v" (w4)
  : [mask] "s" (mask2b), [u1] "v" (u1), [u2] "v" (u2), [u3] "v" (u3), [u4] "v" (u4)
  : "vcc");

  // clang-format on

  /*
    const auto lane = __lane_id() % 4;

    if (lane == 0) {
      const auto v11 = dpp4<0>(v1);
      const auto v12 = dpp4<1>(v1);
      const auto v13 = dpp4<2>(v1);
      const auto v14 = dpp4<3>(v1);
      w1 = v11;
      w2 = v12;
      w3 = v13;
      w4 = v14;
    }
    if (lane == 1) {
      const auto v21 = dpp4<0>(v2);
    const auto v22 = dpp4<1>(v2);
    const auto v23 = dpp4<2>(v2);
    const auto v24 = dpp4<3>(v2);
      w1 = v21;
      w2 = v22;
      w3 = v23;
      w4 = v24;
    }
    if (lane == 2) {
      const auto v31 = dpp4<0>(v3);
    const auto v32 = dpp4<1>(v3);
    const auto v33 = dpp4<2>(v3);
    const auto v34 = dpp4<3>(v3);
      w1 = v31;
      w2 = v32;
      w3 = v33;
      w4 = v34;
    }
    if (lane == 3) {

    const auto v41 = dpp4<0>(v4);
    const auto v42 = dpp4<1>(v4);
    const auto v43 = dpp4<2>(v4);
    const auto v44 = dpp4<3>(v4);
      w1 = v41;
      w2 = v42;
      w3 = v43;
      w4 = v44;
    }*/
}

template <typename T, typename S>
__device__ __forceinline__ auto transpose4x1(T& w, S v) {
  w[0] = dpp4<0>(v);
  w[1] = dpp4<1>(v);
  w[2] = dpp4<2>(v);
  w[3] = dpp4<3>(v);
}

template <typename T>
__device__ __forceinline__ auto transpose4x4(T& w1, T& w2, T& w3, T& w4, T v1, T v2, T v3, T v4) {

  const uint64_t mask1a = 0x5555555555555555ULL;
  const uint64_t mask1b = 0xaaaaaaaaaaaaaaaaULL;
  const uint64_t mask2a = 0x3333333333333333ULL;
  const uint64_t mask2b = 0xccccccccccccccccULL;

  T u1, u2, u3, u4;

  // 11 12 13 14
  // 21 22 23 24
  // 31 32 33 34
  // 41 42 43 44

  // 11 21 13 23 (DPP for row 2)
  // 12 22 14 24 (DPP for row 1)
  // 31 41 33 43 (DPP for row 4)
  // 32 42 34 44 (DPP for row 3)

  // 11 21 31 41 (DPP for row 3)
  // 12 22 32 42 (DPP for row 4)
  // 13 23 33 43 (DPP for row 1)
  // 14 24 34 44 (DPP for row 2)

  // clang-format off

  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CM4STR(0, 0, 2, 2, "%[u1]", "%[v2]", "%[v1]") "\n\t"
  CM4STR(0, 0, 2, 2, "%[u3]", "%[v4]", "%[v3]") "\n\t"
  : [u1] "=v" (u1), [u3] "=v" (u3)
  : [mask] "s" (mask1a), [v1] "v" (v1), [v2] "v" (v2), [v3] "v" (v3), [v4] "v" (v4)
  : "vcc");
  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CM4STR(1, 1, 3, 3, "%[u2]", "%[v1]", "%[v2]") "\n\t"
  CM4STR(1, 1, 3, 3, "%[u4]", "%[v3]", "%[v4]") "\n\t"
  : [u2] "=v" (u2), [u4] "=v" (u4)
  : [mask] "s" (mask1b), [v1] "v" (v1), [v2] "v" (v2), [v3] "v" (v3), [v4] "v" (v4)
  : "vcc");
  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CM4STR(0, 1, 0, 1, "%[w1]", "%[u3]", "%[u1]") "\n\t"
  CM4STR(0, 1, 0, 1, "%[w2]", "%[u4]", "%[u2]") "\n\t"
  : [w1] "=v" (w1), [w2] "=v" (w2)
  : [mask] "s" (mask2a), [u1] "v" (u1), [u2] "v" (u2), [u3] "v" (u3), [u4] "v" (u4)
  : "vcc");
  __asm("s_mov_b64 vcc, %[mask] \n\t"
  CM4STR(2, 3, 2, 3, "%[w3]", "%[u1]", "%[u3]") "\n\t"
  CM4STR(2, 3, 2, 3, "%[w4]", "%[u2]", "%[u4]")
  : [w3] "=v" (w3), [w4] "=v" (w4)
  : [mask] "s" (mask2b), [u1] "v" (u1), [u2] "v" (u2), [u3] "v" (u3), [u4] "v" (u4)
  : "vcc");

  // clang-format on

  /*
    const auto lane = __lane_id() % 4;

    if (lane == 0) {
      const auto v11 = dpp4<0>(v1);
      const auto v12 = dpp4<1>(v1);
      const auto v13 = dpp4<2>(v1);
      const auto v14 = dpp4<3>(v1);
      w1 = v11;
      w2 = v12;
      w3 = v13;
      w4 = v14;
    }
    if (lane == 1) {
      const auto v21 = dpp4<0>(v2);
    const auto v22 = dpp4<1>(v2);
    const auto v23 = dpp4<2>(v2);
    const auto v24 = dpp4<3>(v2);
      w1 = v21;
      w2 = v22;
      w3 = v23;
      w4 = v24;
    }
    if (lane == 2) {
      const auto v31 = dpp4<0>(v3);
    const auto v32 = dpp4<1>(v3);
    const auto v33 = dpp4<2>(v3);
    const auto v34 = dpp4<3>(v3);
      w1 = v31;
      w2 = v32;
      w3 = v33;
      w4 = v34;
    }
    if (lane == 3) {

    const auto v41 = dpp4<0>(v4);
    const auto v42 = dpp4<1>(v4);
    const auto v43 = dpp4<2>(v4);
    const auto v44 = dpp4<3>(v4);
      w1 = v41;
      w2 = v42;
      w3 = v43;
      w4 = v44;
    }*/
}

/*

const auto kdivLocal0 = *(float4*)&kdivCache[k * 256 + threadIdx.x * 4];
const auto kdivLocal1 = *(float4*)&kdivCache[(k + 1) * 256 + threadIdx.x * 4];
const auto kdivLocal2 = *(float4*)&kdivCache[(k + 2) * 256 + threadIdx.x * 4];
const auto kdivLocal3 = *(float4*)&kdivCache[(k + 3) * 256 + threadIdx.x * 4];

float4 kdT[4]{};

transpose4x4(kdT[0].x,
kdT[1].x,
kdT[2].x,
kdT[3].x,
kdivLocal0.x,
kdivLocal1.x,
kdivLocal2.x,
kdivLocal3.x);
transpose4x4(kdT[0].y,
kdT[1].y,
kdT[2].y,
kdT[3].y,
kdivLocal0.y,
kdivLocal1.y,
kdivLocal2.y,
kdivLocal3.y);
transpose4x4(kdT[0].z,
kdT[1].z,
kdT[2].z,
kdT[3].z,
kdivLocal0.z,
kdivLocal1.z,
kdivLocal2.z,
kdivLocal3.z);
transpose4x4(kdT[0].w,
kdT[1].w,
kdT[2].w,
kdT[3].w,
kdivLocal0.w,
kdivLocal1.w,
kdivLocal2.w,
kdivLocal3.w);

#pragma unroll
for (int j = 0; j < Quantities; j += 4) {
  float4 dq4{};
  transpose4x4(dq4.x, dq4.y, dq4.z, dq4.w, dq[j + 0], dq[j + 1], dq[j + 2], dq[j + 3]);

  acc[0][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(kdivLocal.x, dq4.x, acc[0][j / 4], 0, 0, 0);
  acc[1][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(kdivLocal.y, dq4.y, acc[0][j / 4], 0, 0, 0);
  acc[2][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(kdivLocal.z, dq4.z, acc[0][j / 4], 0, 0, 0);
  acc[3][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(kdivLocal.w, dq4.w, acc[0][j / 4], 0, 0, 0);

  #pragma unroll
  for (int kk = 0; kk < 4; ++kk) {
    acc[f][j / 4] =
    __builtin_amdgcn_mfma_f32_4x4x1f32(kdT[jj].x, dq[j + 0], acc[f][j / 4], 0, 0, 0);
    acc[f][j / 4] =
    __builtin_amdgcn_mfma_f32_4x4x1f32(kdT[1].x, dq[j + 1], acc[f][j / 4], 0, 0, 0);
    acc[f][j / 4] =
    __builtin_amdgcn_mfma_f32_4x4x1f32(kdT[2].x, dq[j + 2], acc[f][j / 4], 0, 0, 0);
    acc[f][j / 4] =
    __builtin_amdgcn_mfma_f32_4x4x1f32(kdT[3].x, dq[j + 3], acc[f][j / 4], 0, 0, 0);
}
acc[f][j / 4] =
__builtin_amdgcn_mfma_f32_4x4x1f32(kdT[0].x, dq[j + 0], acc[f][j / 4], 0, 0, 0);
acc[f][j / 4] =
__builtin_amdgcn_mfma_f32_4x4x1f32(kdT[1].x, dq[j + 1], acc[f][j / 4], 0, 0, 0);
acc[f][j / 4] =
__builtin_amdgcn_mfma_f32_4x4x1f32(kdT[2].x, dq[j + 2], acc[f][j / 4], 0, 0, 0);
acc[f][j / 4] =
__builtin_amdgcn_mfma_f32_4x4x1f32(kdT[3].x, dq[j + 3], acc[f][j / 4], 0, 0, 0);
}
*/

template <typename T>
__device__ __forceinline__ auto transpose4x4(af4& w, T v1, T v2, T v3, T v4) {

  // REMARK: we could combine DPP with cndmask (I tried it with inline assembly;
  // but it didn't work w/o errors alas)

  const auto vv2 = dpp<0xa0, 0xf, 0xf, true>(v2);
  const auto vv4 = dpp<0xa0, 0xf, 0xf, true>(v4);
  const auto vv1 = dpp<0xf5, 0xf, 0xf, true>(v1);
  const auto vv3 = dpp<0xf5, 0xf, 0xf, true>(v3);

  const auto u1 = __lane_id() % 2 == 0 ? v1 : vv2;
  const auto u2 = __lane_id() % 2 == 1 ? v2 : vv1;
  const auto u3 = __lane_id() % 2 == 0 ? v3 : vv4;
  const auto u4 = __lane_id() % 2 == 1 ? v4 : vv3;

  const auto uu1 = dpp<0xee, 0xf, 0xf, true>(u1);
  const auto uu2 = dpp<0xee, 0xf, 0xf, true>(u2);
  const auto uu3 = dpp<0x44, 0xf, 0xf, true>(u3);
  const auto uu4 = dpp<0x44, 0xf, 0xf, true>(u4);

  w[0] = __lane_id() % 4 < 2 ? u1 : uu3;
  w[1] = __lane_id() % 4 < 2 ? u2 : uu4;
  w[2] = __lane_id() % 4 >= 2 ? u3 : uu1;
  w[3] = __lane_id() % 4 >= 2 ? u4 : uu2;
}

#define FMADPP4(pos, c, a, b)                                                                      \
  __asm("v_fmac_f32_dpp %0, %1, %2 quad_perm:[" STR(pos) "," STR(pos) "," STR(pos) "," STR(        \
            pos) "] row_mask:0xf bank_mask:0xf bound_ctrl:1"                                       \
        : "+v"(c)                                                                                  \
        : "v"(a), "v"(b)                                                                           \
        :)
#define FMADPP16(pos, c, a, b)                                                                     \
  __asm("v_fmac_f32_dpp %0, %1, %2 row_newbcast:" STR(                                             \
            pos) " row_mask:0xf bank_mask:0xf bound_ctrl:1"                                        \
        : "+v"(c)                                                                                  \
        : "v"(a), "v"(b)                                                                           \
        :)
#define DMADPP16(pos, c, a, b)                                                                     \
  __asm("v_fmac_f64_dpp %0, %1, %2 row_newbcast:" STR(                                             \
            pos) " row_mask:0xf bank_mask:0xf bound_ctrl:1"                                        \
        : "+v"(c)                                                                                  \
        : "v"(a), "v"(b)                                                                           \
        :)

// format:
// c: accumulator
// a: DPP-broadcasted register
// b: multiplicand (vector reg)

template <int Row>
__device__ __forceinline__ void fmacdpp4(float& c, float a, float b);

template <int Row>
__device__ __forceinline__ void fmacdpp16(float& c, float a, float b);

template <int Row>
__device__ __forceinline__ void fmacdpp16(double& c, double a, double b);

template <>
__device__ __forceinline__ void fmacdpp4<0>(float& c, float a, float b) {
  FMADPP4(0, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp4<1>(float& c, float a, float b) {
  FMADPP4(1, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp4<2>(float& c, float a, float b) {
  FMADPP4(2, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp4<3>(float& c, float a, float b) {
  FMADPP4(3, c, a, b);
}

template <>
__device__ __forceinline__ void fmacdpp16<0>(float& c, float a, float b) {
  FMADPP16(0x0, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<1>(float& c, float a, float b) {
  FMADPP16(0x1, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<2>(float& c, float a, float b) {
  FMADPP16(0x2, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<3>(float& c, float a, float b) {
  FMADPP16(0x3, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<4>(float& c, float a, float b) {
  FMADPP16(0x4, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<5>(float& c, float a, float b) {
  FMADPP16(0x5, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<6>(float& c, float a, float b) {
  FMADPP16(0x6, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<7>(float& c, float a, float b) {
  FMADPP16(0x7, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<8>(float& c, float a, float b) {
  FMADPP16(0x8, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<9>(float& c, float a, float b) {
  FMADPP16(0x9, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<10>(float& c, float a, float b) {
  FMADPP16(0xa, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<11>(float& c, float a, float b) {
  FMADPP16(0xb, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<12>(float& c, float a, float b) {
  FMADPP16(0xc, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<13>(float& c, float a, float b) {
  FMADPP16(0xd, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<14>(float& c, float a, float b) {
  FMADPP16(0xe, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<15>(float& c, float a, float b) {
  FMADPP16(0xf, c, a, b);
}

template <int row>
__device__ __forceinline__ void fmacdpp16(double& c, double a, double b);

template <>
__device__ __forceinline__ void fmacdpp16<0>(double& c, double a, double b) {
  DMADPP16(0x0, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<1>(double& c, double a, double b) {
  DMADPP16(0x1, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<2>(double& c, double a, double b) {
  DMADPP16(0x2, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<3>(double& c, double a, double b) {
  DMADPP16(0x3, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<4>(double& c, double a, double b) {
  DMADPP16(0x4, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<5>(double& c, double a, double b) {
  DMADPP16(0x5, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<6>(double& c, double a, double b) {
  DMADPP16(0x6, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<7>(double& c, double a, double b) {
  DMADPP16(0x7, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<8>(double& c, double a, double b) {
  DMADPP16(0x8, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<9>(double& c, double a, double b) {
  DMADPP16(0x9, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<10>(double& c, double a, double b) {
  DMADPP16(0xa, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<11>(double& c, double a, double b) {
  DMADPP16(0xb, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<12>(double& c, double a, double b) {
  DMADPP16(0xc, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<13>(double& c, double a, double b) {
  DMADPP16(0xd, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<14>(double& c, double a, double b) {
  DMADPP16(0xe, c, a, b);
}
template <>
__device__ __forceinline__ void fmacdpp16<15>(double& c, double a, double b) {
  DMADPP16(0xf, c, a, b);
}

template <typename T>
using gptr = __attribute__((address_space(1))) T* __restrict;

template <int k>
__device__ __forceinline__ void local_step(const float* __restrict kdivCache,
                                           float4 dq4[2],
                                           float dq[9],
                                           af4 acc[4][2],
                                           float interm[4][9]) {
  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  const auto kdivLocal = *(const float4* __restrict)&kdivCache[k * 256 + threadIdx.x * 4];

  for (int j = 0; j < (Quantities / 4) * 4; j += 4) {
    acc[0][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.x, acc[0][j / 4], 4, k / 4, 0);
    acc[1][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.y, acc[1][j / 4], 4, k / 4, 0);
    acc[2][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.z, acc[2][j / 4], 4, k / 4, 0);
    acc[3][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.w, acc[3][j / 4], 4, k / 4, 0);
  }

  for (int j = (Quantities / 4) * 4; j < Quantities; ++j) {
    {
      const auto value = readlane(dq[j], k);
      interm[0][j] += kdivLocal.x * value;
      interm[1][j] += kdivLocal.y * value;
      interm[2][j] += kdivLocal.z * value;
      interm[3][j] += kdivLocal.w * value;
    }
  }
}

template <int k>
__device__ __forceinline__ void local_step(
    const float* __restrict kdivCache, af4 dq4[2], float dq[9], af4 acc[4][2], float interm[4][9]) {
  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  const auto kdivLocal = *(const float4* __restrict)&kdivCache[k * 256 + threadIdx.x * 4];

  for (int j = 0; j < (Quantities / 4) * 4; j += 4) {
    acc[0][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.x, acc[0][j / 4], 4, k / 4, 0);
    acc[1][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.y, acc[1][j / 4], 4, k / 4, 0);
    acc[2][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.z, acc[2][j / 4], 4, k / 4, 0);
    acc[3][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.w, acc[3][j / 4], 4, k / 4, 0);
  }

  for (int j = (Quantities / 4) * 4; j < Quantities; ++j) {
    {
      const auto value = readlane(dq[j], k);
      interm[0][j] += kdivLocal.x * value;
      interm[1][j] += kdivLocal.y * value;
      interm[2][j] += kdivLocal.z * value;
      interm[3][j] += kdivLocal.w * value;
    }
  }
}

template <int k>
__device__ __forceinline__ void
    local_step(const float* __restrict kdivCache, af4 dq4[3], af4 acc[4][3]) {
  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  const auto kdivLocal = *(const float4* __restrict)&kdivCache[k * 256 + threadIdx.x * 4];

#pragma unroll
  for (int j = 0; j < 12; j += 4) {
    acc[0][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.x, acc[0][j / 4], 4, k / 4, 0);
    acc[1][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.y, acc[1][j / 4], 4, k / 4, 0);
    acc[2][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.z, acc[2][j / 4], 4, k / 4, 0);
    acc[3][j / 4] = __builtin_amdgcn_mfma_f32_4x4x1f32(
        dq4[j / 4][k % 4], kdivLocal.w, acc[3][j / 4], 4, k / 4, 0);
  }
}

template <int Idx>
__device__ __forceinline__ void hadd(float result[9], float interm[4][9], float star[]) {
  fmacdpp16<Idx % 16>(result[(Idx / 9) % 9], star[Idx / 16], interm[Idx / 81][Idx % 9]);
}

template <int Base, int Idx, int End>
__device__ __forceinline__ void haddMany(float result[9], float interm[4][9], float star[]) {
  if constexpr (Idx < End) {
    hadd<Base + Idx>(result, interm, star);
    haddMany<Base, Idx + 1, End>(result, interm, star);
  }
}

constexpr auto LaunchSize = 256;

__launch_bounds__(LaunchSize) __global__ void kernel_local8(const float** A,
                                                            const float** B,
                                                            unsigned Boffset,
                                                            const float* C1,
                                                            const float* C2,
                                                            const float* C3,
                                                            const float* C4,
                                                            float** D,
                                                            size_t numElements,
                                                            const unsigned* flags) {

  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  __shared__ float kdivCache[64 * 56 * Faces];
  // TODO: maybe try add __shared__ float broadcastCache[64 * 4 * 8]; ?

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 16;
  constexpr int CountR = Count % 16;

  const auto linear = threadIdx.y * blockDim.x + threadIdx.x;

  const float* const __restrict__ glbC1 = C1;
  const float* const __restrict__ glbC2 = C2;
  const float* const __restrict__ glbC3 = C3;
  const float* const __restrict__ glbC4 = C4;

  constexpr int Wavecount = LaunchSize / 64;

#pragma unroll
  for (int i = 0; i < 56 / Wavecount; ++i) {
    auto x1 =
        __builtin_nontemporal_load(&glbC1[i * 56 * Wavecount + threadIdx.y * 56 + threadIdx.x]);
    auto x2 =
        __builtin_nontemporal_load(&glbC2[i * 56 * Wavecount + threadIdx.y * 56 + threadIdx.x]);
    auto x3 =
        __builtin_nontemporal_load(&glbC3[i * 56 * Wavecount + threadIdx.y * 56 + threadIdx.x]);
    auto x4 =
        __builtin_nontemporal_load(&glbC4[i * 56 * Wavecount + threadIdx.y * 56 + threadIdx.x]);

    if (threadIdx.x >= 56) {
      x1 = 0;
      x2 = 0;
      x3 = 0;
      x4 = 0;
    }

    float4 x{};
    x.x = x1;
    x.y = x2;
    x.z = x3;
    x.w = x4;

    *(float4*)&kdivCache[i * 256 * Wavecount + linear * 4] = x;
  }
  __syncthreads();

  for (int b = blockIdx.x * blockDim.y + threadIdx.y; b < numElements;
       b += gridDim.x * blockDim.y) {
    auto glbA = (__attribute__((address_space(1))) const float* const __restrict__)A[b];
    // const float* const __restrict__ glbB = B[b] + Boffset;
    auto glbB = (__attribute__((address_space(1))) const float* const __restrict__)(B[b] + Boffset);
    // const gptr<float> glbB = B[b] + Boffset;
    auto glbD = (__attribute__((address_space(1))) float* const __restrict__)D[b];

    float result[Quantities]{};
    float dq[Quantities]{};
    float star[CountH + 1]{};

    const auto flag = flags[b];

    const bool has1 = (flag & 1) != 0;
    const bool has2 = (flag & 2) != 0;
    const bool has3 = (flag & 4) != 0;
    const bool has4 = (flag & 8) != 0;

    const bool has[4]{has1, has2, has3, has4};

    // load matrices

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      dq[i] = __builtin_nontemporal_load(&glbA[i * 64 + threadIdx.x]);
    }

    float interm[Faces][Quantities]{};

    af4 acc[Faces][3]{};

    af4 dq4[3]{};

    for (int i = 0; i < 2; ++i) {
      transpose4x4(dq4[i], dq[4 * i + 0], dq[4 * i + 1], dq[4 * i + 2], dq[4 * i + 3]);
    }
    transpose4x1(dq4[2], dq[8]);

#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      star[i] = glbB[threadIdx.x % 16 + i * 16];
    }
    if (threadIdx.x % 16 < CountR) {
      star[CountH] = glbB[threadIdx.x % 16 + CountH * 16];
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      result[i] = __builtin_nontemporal_load(&glbD[i * 64 + threadIdx.x]);
    }

    // matmul #1 X = (M @ I) × 4
    local_step<0>(kdivCache, dq4, acc);
    local_step<1>(kdivCache, dq4, acc);
    local_step<2>(kdivCache, dq4, acc);
    local_step<3>(kdivCache, dq4, acc);
    local_step<4>(kdivCache, dq4, acc);
    local_step<5>(kdivCache, dq4, acc);
    local_step<6>(kdivCache, dq4, acc);
    local_step<7>(kdivCache, dq4, acc);
    local_step<8>(kdivCache, dq4, acc);
    local_step<9>(kdivCache, dq4, acc);
    local_step<10>(kdivCache, dq4, acc);
    local_step<11>(kdivCache, dq4, acc);
    local_step<12>(kdivCache, dq4, acc);
    local_step<13>(kdivCache, dq4, acc);
    local_step<14>(kdivCache, dq4, acc);
    local_step<15>(kdivCache, dq4, acc);
    local_step<16>(kdivCache, dq4, acc);
    local_step<17>(kdivCache, dq4, acc);
    local_step<18>(kdivCache, dq4, acc);
    local_step<19>(kdivCache, dq4, acc);
    local_step<20>(kdivCache, dq4, acc);
    local_step<21>(kdivCache, dq4, acc);
    local_step<22>(kdivCache, dq4, acc);
    local_step<23>(kdivCache, dq4, acc);
    local_step<24>(kdivCache, dq4, acc);
    local_step<25>(kdivCache, dq4, acc);
    local_step<26>(kdivCache, dq4, acc);
    local_step<27>(kdivCache, dq4, acc);
    local_step<28>(kdivCache, dq4, acc);
    local_step<29>(kdivCache, dq4, acc);
    local_step<30>(kdivCache, dq4, acc);
    local_step<31>(kdivCache, dq4, acc);
    local_step<32>(kdivCache, dq4, acc);
    local_step<33>(kdivCache, dq4, acc);
    local_step<34>(kdivCache, dq4, acc);
    local_step<35>(kdivCache, dq4, acc);
    local_step<36>(kdivCache, dq4, acc);
    local_step<37>(kdivCache, dq4, acc);
    local_step<38>(kdivCache, dq4, acc);
    local_step<39>(kdivCache, dq4, acc);
    local_step<40>(kdivCache, dq4, acc);
    local_step<41>(kdivCache, dq4, acc);
    local_step<42>(kdivCache, dq4, acc);
    local_step<43>(kdivCache, dq4, acc);
    local_step<44>(kdivCache, dq4, acc);
    local_step<45>(kdivCache, dq4, acc);
    local_step<46>(kdivCache, dq4, acc);
    local_step<47>(kdivCache, dq4, acc);
    local_step<48>(kdivCache, dq4, acc);
    local_step<49>(kdivCache, dq4, acc);
    local_step<50>(kdivCache, dq4, acc);
    local_step<51>(kdivCache, dq4, acc);
    local_step<52>(kdivCache, dq4, acc);
    local_step<53>(kdivCache, dq4, acc);
    local_step<54>(kdivCache, dq4, acc);
    local_step<55>(kdivCache, dq4, acc);

#pragma unroll
    for (int j = 0; j < Quantities; ++j) {
#pragma unroll
      for (int f = 0; f < 4; ++f) {
        interm[f][j] = acc[f][j / 4][j % 4];
      }
    }

    if (has[0]) {
      haddMany<0, 0, 81>(result, interm, star);
    }
    if (has[1]) {
      haddMany<81 * 1, 0, 81>(result, interm, star);
    }
    if (has[2]) {
      haddMany<81 * 2, 0, 81>(result, interm, star);
    }
    if (has[3]) {
      haddMany<81 * 3, 0, 81>(result, interm, star);
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      __builtin_nontemporal_store(result[i], &glbD[i * 64 + threadIdx.x]);
    }
  }
}

__launch_bounds__(LaunchSize) __global__ void kernel_local8b(const float** A,
                                                             const float** B,
                                                             unsigned Boffset,
                                                             const float* C1,
                                                             const float* C2,
                                                             const float* C3,
                                                             const float* C4,
                                                             float** D,
                                                             size_t numElements,
                                                             const unsigned* flags) {

  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  __shared__ float kdivCache[64 * 56 * Faces];
  // TODO: maybe try add __shared__ float broadcastCache[64 * 4 * 8]; ?

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 16;
  constexpr int CountR = Count % 16;

  const auto linear = threadIdx.y * blockDim.x + threadIdx.x;

  const float* const __restrict__ glbC1 = C1;
  const float* const __restrict__ glbC2 = C2;
  const float* const __restrict__ glbC3 = C3;
  const float* const __restrict__ glbC4 = C4;

#pragma unroll
  for (int i = 0; i < 56 / 8; ++i) {
    auto x1 = __builtin_nontemporal_load(&glbC1[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x2 = __builtin_nontemporal_load(&glbC2[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x3 = __builtin_nontemporal_load(&glbC3[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x4 = __builtin_nontemporal_load(&glbC4[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);

    if (threadIdx.x >= 56) {
      x1 = 0;
      x2 = 0;
      x3 = 0;
      x4 = 0;
    }

    float4 x{};
    x.x = x1;
    x.y = x2;
    x.z = x3;
    x.w = x4;

    *(float4*)&kdivCache[i * 256 * 8 + linear * 4] = x;
  }
  __syncthreads();

  for (int b = blockIdx.x * blockDim.y + threadIdx.y; b < numElements;
       b += gridDim.x * blockDim.y) {
    auto glbA = (__attribute__((address_space(1))) const float* const __restrict__)A[b];
    // const float* const __restrict__ glbB = B[b] + Boffset;
    auto glbB = (__attribute__((address_space(1))) const float* const __restrict__)(B[b] + Boffset);
    // const gptr<float> glbB = B[b] + Boffset;
    auto glbD = (__attribute__((address_space(1))) float* const __restrict__)D[b];

    float result[Quantities]{};
    float dq[Quantities]{};
    float star[CountH + 1]{};

    const auto flag = flags[b];

    const bool has1 = (flag & 1) != 0;
    const bool has2 = (flag & 2) != 0;
    const bool has3 = (flag & 4) != 0;
    const bool has4 = (flag & 8) != 0;

    const bool has[4]{has1, has2, has3, has4};

    // load matrices

    /*
#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      dq[i] = __builtin_nontemporal_load(&glbA[i * 64 + threadIdx.x]);
    }
    */

    float interm[Faces][Quantities]{};

    af4 acc[Faces][Quantities / 4]{};

    af4 dq4[2]{};

    /*for (int i = 0; i < 2; ++i) {
      transpose4x4(dq4[i].x,
                   dq4[i].y,
                   dq4[i].z,
                   dq4[i].w,
                   dq[4 * i + 0],
                   dq[4 * i + 1],
                   dq[4 * i + 2],
                   dq[4 * i + 3]);
    }*/
    for (int i = 0; i < 8; i += 4) {
      dq4[i / 4] = __builtin_nontemporal_load((af4*)&glbA[i * 64 + threadIdx.x * 4]);
    }
    for (int i = 8; i < Quantities; ++i) {
      dq[i] = __builtin_nontemporal_load(&glbA[i * 64 + threadIdx.x]);
    }

#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      star[i] = glbB[threadIdx.x % 16 + i * 16];
    }
    if (threadIdx.x % 16 < CountR) {
      star[CountH] = glbB[threadIdx.x % 16 + CountH * 16];
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      result[i] = __builtin_nontemporal_load(&glbD[i * 64 + threadIdx.x]);
    }

    // matmul #1 X = (M @ I) × 4
    local_step<0>(kdivCache, dq4, dq, acc, interm);
    local_step<1>(kdivCache, dq4, dq, acc, interm);
    local_step<2>(kdivCache, dq4, dq, acc, interm);
    local_step<3>(kdivCache, dq4, dq, acc, interm);
    local_step<4>(kdivCache, dq4, dq, acc, interm);
    local_step<5>(kdivCache, dq4, dq, acc, interm);
    local_step<6>(kdivCache, dq4, dq, acc, interm);
    local_step<7>(kdivCache, dq4, dq, acc, interm);
    local_step<8>(kdivCache, dq4, dq, acc, interm);
    local_step<9>(kdivCache, dq4, dq, acc, interm);
    local_step<10>(kdivCache, dq4, dq, acc, interm);
    local_step<11>(kdivCache, dq4, dq, acc, interm);
    local_step<12>(kdivCache, dq4, dq, acc, interm);
    local_step<13>(kdivCache, dq4, dq, acc, interm);
    local_step<14>(kdivCache, dq4, dq, acc, interm);
    local_step<15>(kdivCache, dq4, dq, acc, interm);
    local_step<16>(kdivCache, dq4, dq, acc, interm);
    local_step<17>(kdivCache, dq4, dq, acc, interm);
    local_step<18>(kdivCache, dq4, dq, acc, interm);
    local_step<19>(kdivCache, dq4, dq, acc, interm);
    local_step<20>(kdivCache, dq4, dq, acc, interm);
    local_step<21>(kdivCache, dq4, dq, acc, interm);
    local_step<22>(kdivCache, dq4, dq, acc, interm);
    local_step<23>(kdivCache, dq4, dq, acc, interm);
    local_step<24>(kdivCache, dq4, dq, acc, interm);
    local_step<25>(kdivCache, dq4, dq, acc, interm);
    local_step<26>(kdivCache, dq4, dq, acc, interm);
    local_step<27>(kdivCache, dq4, dq, acc, interm);
    local_step<28>(kdivCache, dq4, dq, acc, interm);
    local_step<29>(kdivCache, dq4, dq, acc, interm);
    local_step<30>(kdivCache, dq4, dq, acc, interm);
    local_step<31>(kdivCache, dq4, dq, acc, interm);
    local_step<32>(kdivCache, dq4, dq, acc, interm);
    local_step<33>(kdivCache, dq4, dq, acc, interm);
    local_step<34>(kdivCache, dq4, dq, acc, interm);
    local_step<35>(kdivCache, dq4, dq, acc, interm);
    local_step<36>(kdivCache, dq4, dq, acc, interm);
    local_step<37>(kdivCache, dq4, dq, acc, interm);
    local_step<38>(kdivCache, dq4, dq, acc, interm);
    local_step<39>(kdivCache, dq4, dq, acc, interm);
    local_step<40>(kdivCache, dq4, dq, acc, interm);
    local_step<41>(kdivCache, dq4, dq, acc, interm);
    local_step<42>(kdivCache, dq4, dq, acc, interm);
    local_step<43>(kdivCache, dq4, dq, acc, interm);
    local_step<44>(kdivCache, dq4, dq, acc, interm);
    local_step<45>(kdivCache, dq4, dq, acc, interm);
    local_step<46>(kdivCache, dq4, dq, acc, interm);
    local_step<47>(kdivCache, dq4, dq, acc, interm);
    local_step<48>(kdivCache, dq4, dq, acc, interm);
    local_step<49>(kdivCache, dq4, dq, acc, interm);
    local_step<50>(kdivCache, dq4, dq, acc, interm);
    local_step<51>(kdivCache, dq4, dq, acc, interm);
    local_step<52>(kdivCache, dq4, dq, acc, interm);
    local_step<53>(kdivCache, dq4, dq, acc, interm);
    local_step<54>(kdivCache, dq4, dq, acc, interm);
    local_step<55>(kdivCache, dq4, dq, acc, interm);

#pragma unroll
    for (int j = 0; j < (Quantities / 4) * 4; ++j) {
#pragma unroll
      for (int f = 0; f < 4; ++f) {
        interm[f][j] = acc[f][j / 4][j % 4];
      }
    }

    if (has[0]) {
      haddMany<0, 0, 81>(result, interm, star);
    }
    if (has[1]) {
      haddMany<81 * 1, 0, 81>(result, interm, star);
    }
    if (has[2]) {
      haddMany<81 * 2, 0, 81>(result, interm, star);
    }
    if (has[3]) {
      haddMany<81 * 3, 0, 81>(result, interm, star);
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      __builtin_nontemporal_store(result[i], &glbD[i * 64 + threadIdx.x]);
    }
  }
}

__launch_bounds__(LaunchSize) __global__ void kernel_local8a(const float** A,
                                                             const float** B,
                                                             unsigned Boffset,
                                                             const float* C1,
                                                             const float* C2,
                                                             const float* C3,
                                                             const float* C4,
                                                             float** D,
                                                             size_t numElements,
                                                             const unsigned* flags) {

  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  __shared__ float kdivCache[64 * 56 * Faces];
  // TODO: maybe try add __shared__ float broadcastCache[64 * 4 * 8]; ?

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 16;
  constexpr int CountR = Count % 16;

  const auto linear = threadIdx.y * blockDim.x + threadIdx.x;

  const float* const __restrict__ glbC1 = C1;
  const float* const __restrict__ glbC2 = C2;
  const float* const __restrict__ glbC3 = C3;
  const float* const __restrict__ glbC4 = C4;

  float dqCache[Quantities]{};

#pragma unroll
  for (int i = 0; i < 56 / 8; ++i) {
    auto x1 = __builtin_nontemporal_load(&glbC1[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x2 = __builtin_nontemporal_load(&glbC2[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x3 = __builtin_nontemporal_load(&glbC3[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x4 = __builtin_nontemporal_load(&glbC4[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);

    if (threadIdx.x >= 56) {
      x1 = 0;
      x2 = 0;
      x3 = 0;
      x4 = 0;
    }

    float4 x{};
    x.x = x1;
    x.y = x2;
    x.z = x3;
    x.w = x4;

    *(float4*)&kdivCache[i * 256 * 8 + linear * 4] = x;
  }
  {
    int b = blockIdx.x * blockDim.y + threadIdx.y;
    if (b < numElements) {
      auto glbA = (__attribute__((address_space(1))) const float* const __restrict__)A[b];
#pragma unroll
      for (int i = 0; i < Quantities; ++i) {
        dqCache[i] = __builtin_nontemporal_load(&glbA[i * 64 + threadIdx.x]);
      }
    }
  }
  __syncthreads();

  for (int b = blockIdx.x * blockDim.y + threadIdx.y; b < numElements;
       b += gridDim.x * blockDim.y) {
    auto glbA = (__attribute__((address_space(1))) const float* const __restrict__)A[b];

    auto bn = b + gridDim.x * blockDim.y;
    if (bn >= numElements) {
      bn = b;
    }
    auto glbAX = (__attribute__((address_space(1))) const float* const __restrict__)A[bn];
    // const float* const __restrict__ glbB = B[b] + Boffset;
    auto glbB = (__attribute__((address_space(1))) const float* const __restrict__)(B[b] + Boffset);
    // const gptr<float> glbB = B[b] + Boffset;
    auto glbD = (__attribute__((address_space(1))) float* const __restrict__)D[b];

    float result[Quantities]{};
    float dq[Quantities]{};
    float star[CountH + 1]{};

    const auto flag = flags[b];

    const bool has1 = (flag & 1) != 0;
    const bool has2 = (flag & 2) != 0;
    const bool has3 = (flag & 4) != 0;
    const bool has4 = (flag & 8) != 0;

    const bool has[4]{has1, has2, has3, has4};

    // load matrices

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      dq[i] = dqCache[i];
    }

    float4 dq4[2]{};

    for (int i = 0; i < 2; ++i) {
      transpose4x4(dq4[i].x,
                   dq4[i].y,
                   dq4[i].z,
                   dq4[i].w,
                   dq[4 * i + 0],
                   dq[4 * i + 1],
                   dq[4 * i + 2],
                   dq[4 * i + 3]);
    }

#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      star[i] = glbB[threadIdx.x % 16 + i * 16];
    }
    if (threadIdx.x % 16 < CountR) {
      star[CountH] = glbB[threadIdx.x % 16 + CountH * 16];
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      result[i] = __builtin_nontemporal_load(&glbD[i * 64 + threadIdx.x]);
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      dqCache[i] = __builtin_nontemporal_load(&glbAX[i * 64 + threadIdx.x]);
    }

    // matmul #1 X = (M @ I) × 4
    float interm[Faces][Quantities]{};

    af4 acc[Faces][Quantities / 4]{};

    local_step<0>(kdivCache, dq4, dq, acc, interm);
    local_step<1>(kdivCache, dq4, dq, acc, interm);
    local_step<2>(kdivCache, dq4, dq, acc, interm);
    local_step<3>(kdivCache, dq4, dq, acc, interm);
    local_step<4>(kdivCache, dq4, dq, acc, interm);
    local_step<5>(kdivCache, dq4, dq, acc, interm);
    local_step<6>(kdivCache, dq4, dq, acc, interm);
    local_step<7>(kdivCache, dq4, dq, acc, interm);
    local_step<8>(kdivCache, dq4, dq, acc, interm);
    local_step<9>(kdivCache, dq4, dq, acc, interm);
    local_step<10>(kdivCache, dq4, dq, acc, interm);
    local_step<11>(kdivCache, dq4, dq, acc, interm);
    local_step<12>(kdivCache, dq4, dq, acc, interm);
    local_step<13>(kdivCache, dq4, dq, acc, interm);
    local_step<14>(kdivCache, dq4, dq, acc, interm);
    local_step<15>(kdivCache, dq4, dq, acc, interm);
    local_step<16>(kdivCache, dq4, dq, acc, interm);
    local_step<17>(kdivCache, dq4, dq, acc, interm);
    local_step<18>(kdivCache, dq4, dq, acc, interm);
    local_step<19>(kdivCache, dq4, dq, acc, interm);
    local_step<20>(kdivCache, dq4, dq, acc, interm);
    local_step<21>(kdivCache, dq4, dq, acc, interm);
    local_step<22>(kdivCache, dq4, dq, acc, interm);
    local_step<23>(kdivCache, dq4, dq, acc, interm);
    local_step<24>(kdivCache, dq4, dq, acc, interm);
    local_step<25>(kdivCache, dq4, dq, acc, interm);
    local_step<26>(kdivCache, dq4, dq, acc, interm);
    local_step<27>(kdivCache, dq4, dq, acc, interm);
    local_step<28>(kdivCache, dq4, dq, acc, interm);
    local_step<29>(kdivCache, dq4, dq, acc, interm);
    local_step<30>(kdivCache, dq4, dq, acc, interm);
    local_step<31>(kdivCache, dq4, dq, acc, interm);
    local_step<32>(kdivCache, dq4, dq, acc, interm);
    local_step<33>(kdivCache, dq4, dq, acc, interm);
    local_step<34>(kdivCache, dq4, dq, acc, interm);
    local_step<35>(kdivCache, dq4, dq, acc, interm);
    local_step<36>(kdivCache, dq4, dq, acc, interm);
    local_step<37>(kdivCache, dq4, dq, acc, interm);
    local_step<38>(kdivCache, dq4, dq, acc, interm);
    local_step<39>(kdivCache, dq4, dq, acc, interm);
    local_step<40>(kdivCache, dq4, dq, acc, interm);
    local_step<41>(kdivCache, dq4, dq, acc, interm);
    local_step<42>(kdivCache, dq4, dq, acc, interm);
    local_step<43>(kdivCache, dq4, dq, acc, interm);
    local_step<44>(kdivCache, dq4, dq, acc, interm);
    local_step<45>(kdivCache, dq4, dq, acc, interm);
    local_step<46>(kdivCache, dq4, dq, acc, interm);
    local_step<47>(kdivCache, dq4, dq, acc, interm);
    local_step<48>(kdivCache, dq4, dq, acc, interm);
    local_step<49>(kdivCache, dq4, dq, acc, interm);
    local_step<50>(kdivCache, dq4, dq, acc, interm);
    local_step<51>(kdivCache, dq4, dq, acc, interm);
    local_step<52>(kdivCache, dq4, dq, acc, interm);
    local_step<53>(kdivCache, dq4, dq, acc, interm);
    local_step<54>(kdivCache, dq4, dq, acc, interm);
    local_step<55>(kdivCache, dq4, dq, acc, interm);

    /*
#pragma unroll //4
    for (int k = 0; k < 56; k += 4) {
      const auto kdivLocal0 = *(const float4* __restrict)&kdivCache[(k + 0) * 256 + threadIdx.x *
4]; const auto kdivLocal1 = *(const float4* __restrict)&kdivCache[(k + 1) * 256 + threadIdx.x * 4];
      const auto kdivLocal2 = *(const float4* __restrict)&kdivCache[(k + 2) * 256 + threadIdx.x *
4]; const auto kdivLocal3 = *(const float4* __restrict)&kdivCache[(k + 3) * 256 + threadIdx.x * 4];

#pragma unroll
      for (int j = 0; j < (Quantities / 4) * 4; j += 4) {
        float4 dq4{};

        dq4.x = transpose4x4bcst(
            dq[j + 0], dq[j + 1], dq[j + 2], dq[j + 3], k + 0, k + 0, k + 0, k + 0);
        dq4.y = transpose4x4bcst(
            dq[j + 0], dq[j + 1], dq[j + 2], dq[j + 3], k + 1, k + 1, k + 1, k + 1);
        dq4.z = transpose4x4bcst(
            dq[j + 0], dq[j + 1], dq[j + 2], dq[j + 3], k + 2, k + 2, k + 2, k + 2);
        dq4.w = transpose4x4bcst(
            dq[j + 0], dq[j + 1], dq[j + 2], dq[j + 3], k + 3, k + 3, k + 3, k + 3);

        acc[0][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.x, kdivLocal0.x, acc[0][j / 4], 0, 0, 0);
        acc[1][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.x, kdivLocal0.y, acc[1][j / 4], 0, 0, 0);
        acc[2][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.x, kdivLocal0.z, acc[2][j / 4], 0, 0, 0);
        acc[3][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.x, kdivLocal0.w, acc[3][j / 4], 0, 0, 0);
        acc[0][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.y, kdivLocal1.x, acc[0][j / 4], 0, 0, 0);
        acc[1][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.y, kdivLocal1.y, acc[1][j / 4], 0, 0, 0);
        acc[2][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.y, kdivLocal1.z, acc[2][j / 4], 0, 0, 0);
        acc[3][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.y, kdivLocal1.w, acc[3][j / 4], 0, 0, 0);
        acc[0][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.z, kdivLocal2.x, acc[0][j / 4], 0, 0, 0);
        acc[1][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.z, kdivLocal2.y, acc[1][j / 4], 0, 0, 0);
        acc[2][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.z, kdivLocal2.z, acc[2][j / 4], 0, 0, 0);
        acc[3][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.z, kdivLocal2.w, acc[3][j / 4], 0, 0, 0);
        acc[0][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.w, kdivLocal3.x, acc[0][j / 4], 0, 0, 0);
        acc[1][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.w, kdivLocal3.y, acc[1][j / 4], 0, 0, 0);
        acc[2][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.w, kdivLocal3.z, acc[2][j / 4], 0, 0, 0);
        acc[3][j / 4] =
            __builtin_amdgcn_mfma_f32_4x4x1f32(dq4.w, kdivLocal3.w, acc[3][j / 4], 0, 0, 0);
      }

#pragma unroll
      for (int j = (Quantities / 4) * 4; j < Quantities; ++j) {
        {
          const auto value = readlane(dq[j], k + 0);
          interm[0][j] += kdivLocal0.x * value;
          interm[1][j] += kdivLocal0.y * value;
          interm[2][j] += kdivLocal0.z * value;
          interm[3][j] += kdivLocal0.w * value;
        }
        {
          const auto value = readlane(dq[j], k + 1);
          interm[0][j] += kdivLocal1.x * value;
          interm[1][j] += kdivLocal1.y * value;
          interm[2][j] += kdivLocal1.z * value;
          interm[3][j] += kdivLocal1.w * value;
        }
        {
          const auto value = readlane(dq[j], k + 2);
          interm[0][j] += kdivLocal2.x * value;
          interm[1][j] += kdivLocal2.y * value;
          interm[2][j] += kdivLocal2.z * value;
          interm[3][j] += kdivLocal2.w * value;
        }
        {
          const auto value = readlane(dq[j], k + 3);
          interm[0][j] += kdivLocal3.x * value;
          interm[1][j] += kdivLocal3.y * value;
          interm[2][j] += kdivLocal3.z * value;
          interm[3][j] += kdivLocal3.w * value;
        }
      }
    }
    */

#pragma unroll
    for (int j = 0; j < (Quantities / 4) * 4; ++j) {
#pragma unroll
      for (int f = 0; f < 4; ++f) {
        interm[f][j] = acc[f][j / 4][j % 4];
        /*transpose4x4(interm[f][j + 0],
                     interm[f][j + 1],
                     interm[f][j + 2],
                     interm[f][j + 3],
                     acc[f][j / 4][0],
                     acc[f][j / 4][1],
                     acc[f][j / 4][2],
                     acc[f][j / 4][3]);*/
      }
    }

    if (has[0]) {
      haddMany<0, 0, 81>(result, interm, star);
    }
    if (has[1]) {
      haddMany<81 * 1, 0, 81>(result, interm, star);
    }
    if (has[2]) {
      haddMany<81 * 2, 0, 81>(result, interm, star);
    }
    if (has[3]) {
      haddMany<81 * 3, 0, 81>(result, interm, star);
    }

    /*
#pragma unroll
    for (int j = 0; j < (Quantities / 4) * 4; ++j) {
#pragma unroll
      for (int f = 0; f < 4; ++f) {
        interm[f][j] = acc[f][j / 4][j % 4];
      }
    }
    */

    // matmul #2 Q += (X @ A*) × 4

    /*
#pragma unroll
    for (int d = 0; d < Faces; ++d) {
      if (has[d]) {
#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
#pragma unroll
          for (int k = 0; k < Quantities; ++k) {
            const auto staridx = k + n * Quantities + Quantities * Quantities * d;
            //result[n] += interm[d][k] * readlane(star[staridx / 64], staridx % 64);
            hadd(result[n], interm[d][k], star[staridx / 16], staridx % 16);
          }
        }
      }
    }*/

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      __builtin_nontemporal_store(result[i], &glbD[i * 64 + threadIdx.x]);
    }
  }
}

__launch_bounds__(LaunchSize) __global__ void kernel_local7b(const float** A,
                                                             const float** B,
                                                             unsigned Boffset,
                                                             const float* C1,
                                                             const float* C2,
                                                             const float* C3,
                                                             const float* C4,
                                                             float** D,
                                                             size_t numElements,
                                                             const unsigned* flags) {

  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  __shared__ float kdivCache[64 * 56 * Faces];
  // TODO: maybe try add __shared__ float broadcastCache[64 * 4 * 8]; ?

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 16;
  constexpr int CountR = Count % 16;

  const auto linear = threadIdx.y * blockDim.x + threadIdx.x;

  const float* const __restrict__ glbC1 = C1;
  const float* const __restrict__ glbC2 = C2;
  const float* const __restrict__ glbC3 = C3;
  const float* const __restrict__ glbC4 = C4;

#pragma unroll
  for (int i = 0; i < 56 / 8; ++i) {
    auto x1 = __builtin_nontemporal_load(&glbC1[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x2 = __builtin_nontemporal_load(&glbC2[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x3 = __builtin_nontemporal_load(&glbC3[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x4 = __builtin_nontemporal_load(&glbC4[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);

    if (threadIdx.x >= 56) {
      x1 = 0;
      x2 = 0;
      x3 = 0;
      x4 = 0;
    }

    float4 x{};
    x.x = x1;
    x.y = x2;
    x.z = x3;
    x.w = x4;

    *(float4*)&kdivCache[i * 256 * 8 + linear * 4] = x;
  }
  __syncthreads();

  for (int b = blockIdx.x * blockDim.y + threadIdx.y; b < numElements;
       b += gridDim.x * blockDim.y) {
    const float* const __restrict__ glbA = A[b];
    const float* const __restrict__ glbB = B[b] + Boffset;
    float* const __restrict__ glbD = D[b];

    float result[Quantities]{};
    float dq[Quantities]{};
    float star[CountH + 1]{};

    const auto flag = flags[b];

    const bool has1 = (flag & 1) != 0;
    const bool has2 = (flag & 2) != 0;
    const bool has3 = (flag & 4) != 0;
    const bool has4 = (flag & 8) != 0;

    const bool has[4]{has1, has2, has3, has4};

    // load matrices

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      dq[i] = __builtin_nontemporal_load(&glbA[i * 64 + threadIdx.x]);
    }

#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      star[i] = glbB[threadIdx.x % 16 + i * 16];
    }
    if (threadIdx.x % 16 < CountR) {
      star[CountH] = glbB[threadIdx.x % 16 + CountH * 16];
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      result[i] = __builtin_nontemporal_load(&glbD[i * 64 + threadIdx.x]);
    }

    // matmul #1 X = (M @ I) × 4
    float interm[Faces][Quantities]{};

#pragma unroll 8
    for (int k = 0; k < 56; ++k) {
      const auto kdivLocal = *(float4*)&kdivCache[k * 256 + threadIdx.x * 4];

#pragma unroll
      for (int j = 0; j < Quantities; ++j) {
        const auto value = readlane(dq[j], k);

        interm[0][j] += kdivLocal.x * value;
        interm[1][j] += kdivLocal.y * value;
        interm[2][j] += kdivLocal.z * value;
        interm[3][j] += kdivLocal.w * value;
      }
    }

    // matmul #2 Q += (X @ A*) × 4
    if (has[0]) {
      haddMany<0, 0, 81>(result, interm, star);
    }
    if (has[1]) {
      haddMany<81 * 1, 0, 81>(result, interm, star);
    }
    if (has[2]) {
      haddMany<81 * 2, 0, 81>(result, interm, star);
    }
    if (has[3]) {
      haddMany<81 * 3, 0, 81>(result, interm, star);
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      __builtin_nontemporal_store(result[i], &glbD[i * 64 + threadIdx.x]);
    }
  }
}

__launch_bounds__(LaunchSize) __global__ void kernel_local7(const float** A,
                                                            const float** B,
                                                            unsigned Boffset,
                                                            const float* C1,
                                                            const float* C2,
                                                            const float* C3,
                                                            const float* C4,
                                                            float** D,
                                                            size_t numElements,
                                                            const unsigned* flags) {

  constexpr int Quantities = 9;
  constexpr int Faces = 4;

  __shared__ float kdivCache[64 * 56 * Faces];
  // TODO: maybe try add __shared__ float broadcastCache[64 * 4 * 8]; ?

  constexpr int Count = Faces * Quantities * Quantities;
  constexpr int CountH = Count / 64;
  constexpr int CountR = Count % 64;

  const auto linear = threadIdx.y * blockDim.x + threadIdx.x;

  const float* const __restrict__ glbC1 = C1;
  const float* const __restrict__ glbC2 = C2;
  const float* const __restrict__ glbC3 = C3;
  const float* const __restrict__ glbC4 = C4;

#pragma unroll
  for (int i = 0; i < 56 / 8; ++i) {
    auto x1 = __builtin_nontemporal_load(&glbC1[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x2 = __builtin_nontemporal_load(&glbC2[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x3 = __builtin_nontemporal_load(&glbC3[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);
    auto x4 = __builtin_nontemporal_load(&glbC4[i * 56 * 8 + threadIdx.y * 56 + threadIdx.x]);

    if (threadIdx.x >= 56) {
      x1 = 0;
      x2 = 0;
      x3 = 0;
      x4 = 0;
    }

    float4 x{};
    x.x = x1;
    x.y = x2;
    x.z = x3;
    x.w = x4;

    *(float4*)&kdivCache[i * 256 * 8 + linear * 4] = x;
  }
  __syncthreads();

  for (int b = blockIdx.x * blockDim.y + threadIdx.y; b < numElements;
       b += gridDim.x * blockDim.y) {
    const float* const __restrict__ glbA = A[b];
    const float* const __restrict__ glbB = B[b] + Boffset;
    float* const __restrict__ glbD = D[b];

    float result[Quantities]{};
    float dq[Quantities]{};
    float star[CountH + 1]{};

    const auto flag = flags[b];

    const bool has1 = (flag & 1) != 0;
    const bool has2 = (flag & 2) != 0;
    const bool has3 = (flag & 4) != 0;
    const bool has4 = (flag & 8) != 0;

    const bool has[4]{has1, has2, has3, has4};

    // load matrices

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      dq[i] = __builtin_nontemporal_load(&glbA[i * 64 + threadIdx.x]);
    }

#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      star[i] = glbB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < CountR) {
      star[CountH] = glbB[threadIdx.x + CountH * 64];
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      result[i] = __builtin_nontemporal_load(&glbD[i * 64 + threadIdx.x]);
    }

    // matmul #1 X = (M @ I) × 4
    float interm[Faces][Quantities]{};

#pragma unroll 8
    for (int k = 0; k < 56; ++k) {
      const auto kdivLocal = *(float4*)&kdivCache[k * 256 + threadIdx.x * 4];

#pragma unroll
      for (int j = 0; j < Quantities; ++j) {
        const auto value = readlane(dq[j], k);

        interm[0][j] += kdivLocal.x * value;
        interm[1][j] += kdivLocal.y * value;
        interm[2][j] += kdivLocal.z * value;
        interm[3][j] += kdivLocal.w * value;
      }
    }

    // matmul #2 Q += (X @ A*) × 4
#pragma unroll
    for (int d = 0; d < Faces; ++d) {
      if (has[d]) {
#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
#pragma unroll
          for (int k = 0; k < Quantities; ++k) {
            const auto staridx = k + n * Quantities + Quantities * Quantities * d;
            result[n] += interm[d][k] * readlane(star[staridx / 64], staridx % 64);
          }
        }
      }
    }

#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      __builtin_nontemporal_store(result[i], &glbD[i * 64 + threadIdx.x]);
    }
  }
}

__launch_bounds__(256) __global__ void kernel_local6(const float** A,
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

  __shared__ __align__(8) float cache1[8 * 56 * 4 + 16];

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;
  const auto linear_x = threadIdx.x + blockDim.x * threadIdx.y;

  // TODO: handle overhang
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

    float reg0[Faces][Quantities][2]{};
    float reg1[Quantities][2]{};

    float* shrmem0 = &total_shrmem0[(576 + Count) * threadIdx.y];

    // writing to shr mem: from reg0 to _1
    float* __restrict__ _1 = &shrmem0[0];
#pragma unroll
    for (int i = 0; i < Quantities; ++i) {
      _1[tid_x + i * 64] = __builtin_nontemporal_load(&glbA[tid_x + i * 64]);
      _1[tid_x + i * 64 + 32] = __builtin_nontemporal_load(&glbA[tid_x + i * 64 + 32]);
    }

    float* __restrict__ _0 = &shrmem0[576];

    // load all 4 9×9 matrices
    // loading glbB to _0: # no trans, extended
#pragma unroll
    for (int i = 0; i < CountH; ++i) {
      _0[tid_x + i * 32] = glbB[tid_x + i * 32];
    }
    if (tid_x < CountR) {
      _0[tid_x + CountH * 32] = glbB[tid_x + CountH * 32];
    }

#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      reg1[n][0] = __builtin_nontemporal_load(&glbD[tid_x + n * 64 + 0]);
      reg1[n][1] = __builtin_nontemporal_load(&glbD[tid_x + n * 64 + 32]);
    }

    __builtin_amdgcn_sched_barrier(0);

    constexpr int Cache = 8;

    // gemm: glbC x _1

    // #pragma unroll
    for (int k = 0; k < 56; k += Cache) {
      if (linear_x * 2 < 448) {
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 0] = *(float2*)&glbC1[linear_x * 2 + k * 56];
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 1] = *(float2*)&glbC2[linear_x * 2 + k * 56];
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 2] = *(float2*)&glbC3[linear_x * 2 + k * 56];
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 3] = *(float2*)&glbC4[linear_x * 2 + k * 56];
        __builtin_amdgcn_sched_group_barrier(0x020, 4, 1);
        __builtin_amdgcn_sched_group_barrier(0x0200, 4, 1);
      }
      __syncthreads();

      float values[Faces][Cache][2]{};
      if (has1) {
#pragma unroll
        for (int kk = 0; kk < Cache; ++kk) {
          values[0][kk][0] = cache1[tid_x + 56 * kk + Cache * 56 * 0];
          values[0][kk][1] = cache1[tid_x + 56 * kk + Cache * 56 * 0 + 32];
        }
        __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
      }
      if (has2) {
#pragma unroll
        for (int kk = 0; kk < Cache; ++kk) {
          values[1][kk][0] = cache1[tid_x + 56 * kk + Cache * 56 * 1];
          values[1][kk][1] = cache1[tid_x + 56 * kk + Cache * 56 * 1 + 32];
        }
        __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
      }
      if (has3) {
#pragma unroll
        for (int kk = 0; kk < Cache; ++kk) {
          values[2][kk][0] = cache1[tid_x + 56 * kk + Cache * 56 * 2];
          values[2][kk][1] = cache1[tid_x + 56 * kk + Cache * 56 * 2 + 32];
        }
        __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
      }
      if (has4) {
#pragma unroll
        for (int kk = 0; kk < Cache; ++kk) {
          values[3][kk][0] = cache1[tid_x + 56 * kk + Cache * 56 * 3];
          values[3][kk][1] = cache1[tid_x + 56 * kk + Cache * 56 * 3 + 32];
        }
        __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
      }

#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
        float local[Cache]{};
#pragma unroll
        for (int kk = 0; kk < Cache; kk += 4) {
          const auto prelocal = *(float4*)&_1[(k + kk) + n * 64];
          local[kk] = prelocal.x;
          local[kk + 1] = prelocal.y;
          local[kk + 2] = prelocal.z;
          local[kk + 3] = prelocal.w;
        }
#pragma unroll
        for (int kk = 0; kk < Cache; ++kk) {
#pragma unroll
          for (int d = 0; d < Faces; ++d) {
            reg0[d][n][0] += values[d][kk][0] * local[kk];
            reg0[d][n][1] += values[d][kk][1] * local[kk];
          }
        }

        __builtin_amdgcn_sched_group_barrier(0x0100, Cache / 4, 0);
      }
    }

#pragma unroll
    for (int x = 0; x < Faces * Quantities * Quantities; x += 4) {
      float4 _0L = *((float4*)&_0[x]);
#pragma unroll
      for (int z = 0; z < 2; ++z) {
        {
          const auto d = x / (Quantities * Quantities);
          const auto n = (x / Quantities) % Quantities;
          const auto k = x % Quantities;
          reg1[n][z] += reg0[d][k][z] * _0L.x;
        }
        {
          const auto d = (x + 1) / (Quantities * Quantities);
          const auto n = ((x + 1) / Quantities) % Quantities;
          const auto k = (x + 1) % Quantities;
          reg1[n][z] += reg0[d][k][z] * _0L.y;
        }
        {
          const auto d = (x + 2) / (Quantities * Quantities);
          const auto n = ((x + 2) / Quantities) % Quantities;
          const auto k = (x + 2) % Quantities;
          reg1[n][z] += reg0[d][k][z] * _0L.z;
        }
        {
          const auto d = (x + 3) / (Quantities * Quantities);
          const auto n = ((x + 3) / Quantities) % Quantities;
          const auto k = (x + 3) % Quantities;
          reg1[n][z] += reg0[d][k][z] * _0L.w;
        }
      }
      __builtin_amdgcn_sched_group_barrier(0x0100, 1, 0);
    }

    // write results back to glb. memory
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg1[n][0], &glbD[tid_x + n * 64]);
      __builtin_nontemporal_store(reg1[n][1], &glbD[tid_x + n * 64 + 32]);
    }
  }
}

__launch_bounds__(256) __global__ void kernel_local5(const float** A,
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

  __shared__ __align__(8) float total_shrmem0[(576 + Count) * 4];

  __shared__ __align__(8) float cache1[8 * 56 * 4];

  const auto tid_x = threadIdx.x;
  unsigned batchId = threadIdx.y + blockDim.y * blockIdx.x;
  const auto linear_x = threadIdx.x + blockDim.x * threadIdx.y;

  // TODO: handle overhang
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

    __builtin_amdgcn_sched_group_barrier(0x220, 2 * Quantities + CountH, 0);

    constexpr int Cache = 8;

    // gemm: glbC x _1

    // #pragma unroll
    for (int k = 0; k < 56; k += Cache) {
      if (linear_x * 2 < 448) {
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 0] = *(float2*)&glbC1[linear_x * 2 + k * 56];
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 1] = *(float2*)&glbC2[linear_x * 2 + k * 56];
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 2] = *(float2*)&glbC3[linear_x * 2 + k * 56];
        *(float2*)&cache1[linear_x * 2 + Cache * 56 * 3] = *(float2*)&glbC4[linear_x * 2 + k * 56];
        __builtin_amdgcn_sched_group_barrier(0x020, 4, 1);
        __builtin_amdgcn_sched_group_barrier(0x0200, 4, 1);
      }
      __syncthreads();

      if (tid_x < 56) {
        float values[Faces][Cache]{};
        if (has1) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[0][kk] = cache1[tid_x + 56 * kk + Cache * 56 * 0];
          }
          __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
        }
        if (has2) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[1][kk] = cache1[tid_x + 56 * kk + Cache * 56 * 1];
          }
          __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
        }
        if (has3) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[2][kk] = cache1[tid_x + 56 * kk + Cache * 56 * 2];
          }
          __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
        }
        if (has4) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[3][kk] = cache1[tid_x + 56 * kk + Cache * 56 * 3];
          }
          __builtin_amdgcn_sched_group_barrier(0x0100, Cache, 1);
        }

#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
          float local[Cache]{};
#pragma unroll
          for (int kk = 0; kk < Cache; kk += 4) {
            const auto prelocal = *(float4*)&_1[(k + kk) + n * 64];
            local[kk] = prelocal.x;
            local[kk + 1] = prelocal.y;
            local[kk + 2] = prelocal.z;
            local[kk + 3] = prelocal.w;
          }
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
#pragma unroll
            for (int d = 0; d < Faces; ++d) {
              reg0[d][n] += values[d][kk] * local[kk];
            }
          }

          __builtin_amdgcn_sched_group_barrier(0x0100, Cache / 4, 0);
        }
      }
    }

#pragma unroll
    for (int x = 0; x < Faces * Quantities * Quantities; x += 4) {
      float4 _0L = *((float4*)&_0[x]);
      {
        const auto d = x / (Quantities * Quantities);
        const auto n = (x / Quantities) % Quantities;
        const auto k = x % Quantities;
        reg1[n] += reg0[d][k] * _0L.x;
      }
      {
        const auto d = (x + 1) / (Quantities * Quantities);
        const auto n = ((x + 1) / Quantities) % Quantities;
        const auto k = (x + 1) % Quantities;
        reg1[n] += reg0[d][k] * _0L.y;
      }
      {
        const auto d = (x + 2) / (Quantities * Quantities);
        const auto n = ((x + 2) / Quantities) % Quantities;
        const auto k = (x + 2) % Quantities;
        reg1[n] += reg0[d][k] * _0L.z;
      }
      {
        const auto d = (x + 3) / (Quantities * Quantities);
        const auto n = ((x + 3) / Quantities) % Quantities;
        const auto k = (x + 3) % Quantities;
        reg1[n] += reg0[d][k] * _0L.w;
      }
      __builtin_amdgcn_sched_group_barrier(0x0100, 1, 0);
    }

    // write results back to glb. memory
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg1[n], &glbD[tid_x + n * 64]);
    }
  }
}

__launch_bounds__(64) __global__ void kernel_local4(const float** A,
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

  __shared__ __align__(8) float total_shrmem0[(576 + Count) * 1];

  __shared__ __align__(8) float cache1[32 * 64];
  __shared__ __align__(8) float cache2[32 * 64];

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

    //__builtin_amdgcn_sched_group_barrier(0x220, Quantities + CountH + (CountR > 0), 0);

    constexpr int Cache = 8;

    __builtin_amdgcn_sched_barrier(0);

    int k = 0;
    if (has1) {
#pragma unroll
      for (int kk = 0; kk < Cache; kk += 4) {
        *(float4*)&cache1[tid_x + 64 * (Cache * 0 + kk)] = *(float4*)&glbC1[tid_x + (k + kk) * 56];
      }
      //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 1);
    }
    if (has2) {
#pragma unroll
      for (int kk = 0; kk < Cache; kk += 4) {
        *(float4*)&cache1[tid_x + 64 * (Cache * 1 + kk)] = *(float4*)&glbC2[tid_x + (k + kk) * 56];
      }
      //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 2);
    }
    if (has3) {
#pragma unroll
      for (int kk = 0; kk < Cache; kk += 4) {
        *(float4*)&cache1[tid_x + 64 * (Cache * 2 + kk)] = *(float4*)&glbC3[tid_x + (k + kk) * 56];
      }
      //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 3);
    }
    if (has4) {
#pragma unroll
      for (int kk = 0; kk < Cache; kk += 4) {
        *(float4*)&cache1[tid_x + 64 * (Cache * 3 + kk)] = *(float4*)&glbC4[tid_x + (k + kk) * 56];
      }
      //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 4);
    }

    __builtin_amdgcn_sched_barrier(0x20);
    __builtin_amdgcn_sched_barrier(0x200);

    float* cache = cache1;
    float* cacheOther = cache2;

    // gemm: glbC x _1
    if (tid_x < 56) {

      // #pragma unroll
      for (int k = 0; k < 56; k += Cache) {
        if (k + Cache < 56) {
          if (has1) {
#pragma unroll
            for (int kk = 0; kk < Cache; kk += 4) {
              *(float4*)&cacheOther[tid_x + 64 * (Cache * 0 + kk)] =
                  *(float4*)&glbC1[tid_x + (k + kk) * 56];
            }
            //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 1);
          }
          if (has2) {
#pragma unroll
            for (int kk = 0; kk < Cache; kk += 4) {
              *(float4*)&cacheOther[tid_x + 64 * (Cache * 1 + kk)] =
                  *(float4*)&glbC2[tid_x + (k + kk) * 56];
            }
            //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 2);
          }
          if (has3) {
#pragma unroll
            for (int kk = 0; kk < Cache; kk += 4) {
              *(float4*)&cacheOther[tid_x + 64 * (Cache * 2 + kk)] =
                  *(float4*)&glbC3[tid_x + (k + kk) * 56];
            }
            //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 3);
          }
          if (has4) {
#pragma unroll
            for (int kk = 0; kk < Cache; kk += 4) {
              *(float4*)&cacheOther[tid_x + 64 * (Cache * 3 + kk)] =
                  *(float4*)&glbC4[tid_x + (k + kk) * 56];
            }
            //__builtin_amdgcn_sched_group_barrier(0x0020, Cache, 4);
          }
        }

        float values[Cache][Faces]{};
#pragma unroll
        for (int kk = 0; kk < Cache; ++kk) {
#pragma unroll
          for (int d = 0; d < Faces; ++d) {
            values[kk][d] = cache1[tid_x + 64 * (Cache * d + kk)];
          }
        }

        __builtin_amdgcn_sched_group_barrier(0x0100, Cache * Faces, 0);

#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
          float local[Cache]{};
#pragma unroll
          for (int kk = 0; kk < Cache; kk += 4) {
            const auto prelocal = *(float4*)&_1[(k + kk) + n * 64];
            local[kk] = prelocal.x;
            local[kk + 1] = prelocal.y;
            local[kk + 2] = prelocal.z;
            local[kk + 3] = prelocal.w;
          }
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
#pragma unroll
            for (int d = 0; d < Faces; ++d) {
              reg0[d][n] += values[kk][d] * local[kk];
            }
          }

          //__builtin_amdgcn_sched_group_barrier(0x2, Cache * Faces, 0);
          __builtin_amdgcn_sched_group_barrier(0x0100, Cache / 4, 0);
        }

        __builtin_amdgcn_sched_barrier(0x20);
        __builtin_amdgcn_sched_barrier(0x200);
        /*
        #pragma unroll
                for (int f = 0; f < Faces; ++f) {
        #pragma unroll
                  for (int kk = 0; kk < Cache; ++kk) {
                    cache1[tid_x + 64*(Cache * f + kk)] = cache2[tid_x + 64*(Cache * f + kk)];
                  }
                }

                __builtin_amdgcn_sched_group_barrier(0x0100, Faces * Cache, 0);
                __builtin_amdgcn_sched_group_barrier(0x0200, Faces * Cache, 0);*/

        float* third = cache;
        cache = cacheOther;
        cacheOther = third;

        /*__builtin_amdgcn_sched_barrier(0x2);*/

        /*if (has1) {
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 0);
        }
        if (has2) {
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 0);
        }
        if (has3) {
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 0);
        }
        if (has4) {
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 0);
        }*/
      }
    }

    __builtin_amdgcn_sched_barrier(0);

#pragma unroll
    for (int x = 0; x < Faces * Quantities * Quantities; x += 4) {
      float4 _0L = *((float4*)&_0[x]);
      {
        const auto d = x / (Quantities * Quantities);
        const auto n = (x / Quantities) % Quantities;
        const auto k = x % Quantities;
        reg1[n] += reg0[d][k] * _0L.x;
      }
      {
        const auto d = (x + 1) / (Quantities * Quantities);
        const auto n = ((x + 1) / Quantities) % Quantities;
        const auto k = (x + 1) % Quantities;
        reg1[n] += reg0[d][k] * _0L.y;
      }
      {
        const auto d = (x + 2) / (Quantities * Quantities);
        const auto n = ((x + 2) / Quantities) % Quantities;
        const auto k = (x + 2) % Quantities;
        reg1[n] += reg0[d][k] * _0L.z;
      }
      {
        const auto d = (x + 3) / (Quantities * Quantities);
        const auto n = ((x + 3) / Quantities) % Quantities;
        const auto k = (x + 3) % Quantities;
        reg1[n] += reg0[d][k] * _0L.w;
      }
      __builtin_amdgcn_sched_group_barrier(0x0100, 1, 0);
    }

    // write results back to glb. memory
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg1[n], &glbD[tid_x + n * 64]);
    }
  }
}

__launch_bounds__(64) __global__ void kernel_local3(const float** A,
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

  __shared__ __align__(8) float total_shrmem0[(576 + Count) * 1];

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

    __builtin_amdgcn_sched_group_barrier(0x220, Quantities + CountH + (CountR > 0), 0);

    constexpr int Cache = 8;

    // gemm: glbC x _1
    if (tid_x < 56) {

      // #pragma unroll
      for (int k = 0; k < 56; k += Cache) {
        float values[Faces][Cache]{};
        if (has1) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[0][kk] = glbC1[tid_x + (k + kk) * 56];
          }
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 1);
        }
        if (has2) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[1][kk] = glbC2[tid_x + (k + kk) * 56];
          }
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 1);
        }
        if (has3) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[2][kk] = glbC3[tid_x + (k + kk) * 56];
          }
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 1);
        }
        if (has4) {
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
            values[3][kk] = glbC4[tid_x + (k + kk) * 56];
          }
          __builtin_amdgcn_sched_group_barrier(0x0020, Cache, 1);
        }

#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
          float local[Cache]{};
#pragma unroll
          for (int kk = 0; kk < Cache; kk += 4) {
            const auto prelocal = *(float4*)&_1[(k + kk) + n * 64];
            local[kk] = prelocal.x;
            local[kk + 1] = prelocal.y;
            local[kk + 2] = prelocal.z;
            local[kk + 3] = prelocal.w;
          }
#pragma unroll
          for (int kk = 0; kk < Cache; ++kk) {
#pragma unroll
            for (int d = 0; d < Faces; ++d) {
              reg0[d][n] += values[d][kk] * local[kk];
            }
          }

          __builtin_amdgcn_sched_group_barrier(0x0100, Cache / 4, 0);
        }
      }
    }

    /*
    // gemm: glbA x _0
#pragma unroll
    for (int d = 0; d < Faces; ++d) {
#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
#pragma unroll
        for (int k = 0; k < Quantities; ++k) {
          reg1[n] += reg0[d][k] * _0[k + n * Quantities + Quantities * Quantities * d];
        }
        __builtin_amdgcn_sched_group_barrier(0x0100, 9, 0);
      }
      // __builtin_amdgcn_sched_group_barrier(0x0100, 8, 0);
    }
      */
#pragma unroll
    for (int x = 0; x < Faces * Quantities * Quantities; x += 4) {
      float4 _0L = *((float4*)&_0[x]);
      {
        const auto d = x / (Quantities * Quantities);
        const auto n = (x / Quantities) % Quantities;
        const auto k = x % Quantities;
        reg1[n] += reg0[d][k] * _0L.x;
      }
      {
        const auto d = (x + 1) / (Quantities * Quantities);
        const auto n = ((x + 1) / Quantities) % Quantities;
        const auto k = (x + 1) % Quantities;
        reg1[n] += reg0[d][k] * _0L.y;
      }
      {
        const auto d = (x + 2) / (Quantities * Quantities);
        const auto n = ((x + 2) / Quantities) % Quantities;
        const auto k = (x + 2) % Quantities;
        reg1[n] += reg0[d][k] * _0L.z;
      }
      {
        const auto d = (x + 3) / (Quantities * Quantities);
        const auto n = ((x + 3) / Quantities) % Quantities;
        const auto k = (x + 3) % Quantities;
        reg1[n] += reg0[d][k] * _0L.w;
      }
      __builtin_amdgcn_sched_group_barrier(0x0100, 1, 0);
    }

    // write results back to glb. memory
#pragma unroll
    for (int n = 0; n < Quantities; ++n) {
      __builtin_nontemporal_store(reg1[n], &glbD[tid_x + n * 64]);
    }
  }
}

__launch_bounds__(512) __global__ void kernel_local2(const float** A,
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

    const float* __restrict__ ptrC[4]{};
    int ptrMap[4]{};

    int ptrcnt = 0;
    if (has1) {
      ptrC[ptrcnt] = glbC1;
      ptrMap[ptrcnt] = 0;
      ++ptrcnt;
    }
    if (has2) {
      ptrC[ptrcnt] = glbC2;
      ptrMap[ptrcnt] = 1;
      ++ptrcnt;
    }
    if (has3) {
      ptrC[ptrcnt] = glbC3;
      ptrMap[ptrcnt] = 2;
      ++ptrcnt;
    }
    if (has4) {
      ptrC[ptrcnt] = glbC4;
      ptrMap[ptrcnt] = 3;
      ++ptrcnt;
    }

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
        if (ptrcnt >= 1) {
#pragma unroll
          for (int kk = 0; kk < 8; ++kk) {
            values[0][kk] = ptrC[0][tid_x + (k + kk) * 56];
          }
          if (ptrcnt >= 2) {
#pragma unroll
            for (int kk = 0; kk < 8; ++kk) {
              values[1][kk] = ptrC[1][tid_x + (k + kk) * 56];
            }
            if (ptrcnt >= 3) {
#pragma unroll
              for (int kk = 0; kk < 8; ++kk) {
                values[2][kk] = ptrC[2][tid_x + (k + kk) * 56];
              }
              if (ptrcnt == 4) {
#pragma unroll
                for (int kk = 0; kk < 8; ++kk) {
                  values[3][kk] = ptrC[3][tid_x + (k + kk) * 56];
                }
              }
            }
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
            if (ptrcnt >= 1) {
              reg0[0][n] += values[0][kk] * local[kk];
              if (ptrcnt >= 2) {
                reg0[1][n] += values[1][kk] * local[kk];
                if (ptrcnt >= 3) {
                  reg0[2][n] += values[2][kk] * local[kk];
                  if (ptrcnt >= 4) {
                    reg0[3][n] += values[3][kk] * local[kk];
                  }
                }
              }
            }
          }
        }
      }
    }

    // gemm: glbA x _0
    if (ptrcnt >= 1) {
#pragma unroll
      for (int n = 0; n < Quantities; ++n) {
#pragma unroll
        for (int k = 0; k < Quantities; ++k) {
          reg1[n] += reg0[0][k] * _0[k + n * Quantities + Quantities * Quantities * ptrMap[0]];
        }
      }
      if (ptrcnt >= 2) {
#pragma unroll
        for (int n = 0; n < Quantities; ++n) {
#pragma unroll
          for (int k = 0; k < Quantities; ++k) {
            reg1[n] += reg0[1][k] * _0[k + n * Quantities + Quantities * Quantities * ptrMap[1]];
          }
        }
        if (ptrcnt >= 3) {
#pragma unroll
          for (int n = 0; n < Quantities; ++n) {
#pragma unroll
            for (int k = 0; k < Quantities; ++k) {
              reg1[n] += reg0[2][k] * _0[k + n * Quantities + Quantities * Quantities * ptrMap[2]];
            }
          }
          if (ptrcnt >= 4) {
#pragma unroll
            for (int n = 0; n < Quantities; ++n) {
#pragma unroll
              for (int k = 0; k < Quantities; ++k) {
                reg1[n] +=
                    reg0[3][k] * _0[k + n * Quantities + Quantities * Quantities * ptrMap[3]];
              }
            }
          }
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
  /*
  dim3 block(64, 1, 1);
  dim3 grid((numElements + 1 - 1) / 1, 1, 1);
  hipStream_t stream = (streamPtr != nullptr) ? static_cast<hipStream_t>(streamPtr) : 0;
  kernel_local3<<<grid, block, 0, stream>>>(A, B, Boffset, C1, C2, C3, C4, D, numElements, flags);
  */

  static int gridsize = -1;

  if (gridsize < 0) {
    int device{}, smCount{}, blocksPerSM{};
    hipGetDevice(&device);
    hipDeviceGetAttribute(&smCount, hipDeviceAttributeMultiprocessorCount, device);
    hipOccupancyMaxActiveBlocksPerMultiprocessor(&blocksPerSM, kernel_local7, LaunchSize, 0);
    if (blocksPerSM > 0) {
      gridsize = smCount * blocksPerSM;
    } else {
      gridsize = smCount;
    }
  }

  dim3 block(64, LaunchSize / 64, 1);
  dim3 grid(gridsize, 1, 1);
  hipStream_t stream = (streamPtr != nullptr) ? static_cast<hipStream_t>(streamPtr) : 0;
  kernel_local7<<<grid, block, 0, stream>>>(A, B, Boffset, C1, C2, C3, C4, D, numElements, flags);
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

  if constexpr (Size1 % 8 != 0) {
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

    if constexpr (Size1 % 8 != 0) {
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
