// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/IntegralAux.h"

#include <algorithm>
#include <cstddef>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

namespace seissol::kernels::device::aux::integral {

namespace {
constexpr unsigned Threads = 128;
constexpr std::size_t MaxBlocks = 65535;

// one block of threads per block of values
template <typename T>
__global__ void __launch_bounds__(Threads) accumulateKernel(const T* const* source,
                                                            double* const* target,
                                                            std::size_t size,
                                                            std::size_t count) {
  for (std::size_t block = blockIdx.x; block < count; block += gridDim.x) {
    const T* __restrict from = source[block];
    double* __restrict to = target[block];
    for (std::size_t i = threadIdx.x; i < size; i += blockDim.x) {
      to[i] += static_cast<double>(from[i]);
    }
  }
}

template <typename T>
void launch(const T* const* source,
            double* const* target,
            std::size_t size,
            std::size_t count,
            void* stream) {
  if (count == 0 || size == 0) {
    return;
  }
  const auto blocks = std::min(count, MaxBlocks);
  accumulateKernel<T>
      <<<blocks, Threads, 0, reinterpret_cast<StreamT>(stream)>>>(source, target, size, count);
}
} // namespace

void accumulate(const float* const* source,
                double* const* target,
                std::size_t size,
                std::size_t count,
                void* stream) {
  launch(source, target, size, count, stream);
}

void accumulate(const double* const* source,
                double* const* target,
                std::size_t size,
                std::size_t count,
                void* stream) {
  launch(source, target, size, count, stream);
}

} // namespace seissol::kernels::device::aux::integral
