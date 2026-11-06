// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/HostCopyAux.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

namespace seissol::kernels::device::aux::hostcopy {

namespace {
constexpr unsigned Threads = 256;
constexpr std::size_t MaxBlocks = 1024;

// in words of 16 bytes where both sides are aligned to them, and the rest in bytes
__global__ void copyToHostKernel(uint4* __restrict target,
                                 const uint4* __restrict source,
                                 std::size_t words,
                                 unsigned char* __restrict targetRest,
                                 const unsigned char* __restrict sourceRest,
                                 std::size_t rest) {
  const std::size_t stride = static_cast<std::size_t>(gridDim.x) * blockDim.x;
  const std::size_t first = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  for (std::size_t i = first; i < words; i += stride) {
    target[i] = source[i];
  }
  for (std::size_t i = first; i < rest; i += stride) {
    targetRest[i] = sourceRest[i];
  }
}
} // namespace

void copyToHost(void* hostOnDevice, const void* device, std::size_t bytes, void* stream) {
  if (bytes == 0) {
    return;
  }
  const bool aligned =
      (reinterpret_cast<std::uintptr_t>(hostOnDevice) | reinterpret_cast<std::uintptr_t>(device)) %
          sizeof(uint4) ==
      0;
  const std::size_t words = aligned ? bytes / sizeof(uint4) : 0;
  const std::size_t rest = bytes - words * sizeof(uint4);
  const auto blocks = std::clamp<std::size_t>(
      (std::max(words, rest) + Threads - 1) / Threads, std::size_t{1}, MaxBlocks);
  copyToHostKernel<<<blocks, Threads, 0, reinterpret_cast<StreamT>(stream)>>>(
      static_cast<uint4*>(hostOnDevice),
      static_cast<const uint4*>(device),
      words,
      static_cast<unsigned char*>(hostOnDevice) + words * sizeof(uint4),
      static_cast<const unsigned char*>(device) + words * sizeof(uint4),
      rest);
}

} // namespace seissol::kernels::device::aux::hostcopy
