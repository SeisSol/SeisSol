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
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::hostcopy {

void copyToHost(void* hostOnDevice, const void* device, std::size_t bytes, void* stream) {
  if (bytes == 0) {
    return;
  }
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  // as the CUDA/HIP kernel: in words of 16 bytes where both sides are aligned to them, and the rest
  // in bytes. A work-item per byte wrote the host memory a byte at a time (on a Max 1550, at about
  // 0.2 GB/s instead of 50).
  using Word = sycl::vec<std::uint32_t, 4>;
  const bool aligned =
      ((reinterpret_cast<std::uintptr_t>(hostOnDevice) | reinterpret_cast<std::uintptr_t>(device)) %
       sizeof(Word)) == 0;
  const std::size_t words = aligned ? bytes / sizeof(Word) : 0;
  const std::size_t rest = bytes - words * sizeof(Word);
  auto* target = static_cast<Word*>(hostOnDevice);
  const auto* source = static_cast<const Word*>(device);
  auto* targetRest = static_cast<unsigned char*>(hostOnDevice) + words * sizeof(Word);
  const auto* sourceRest = static_cast<const unsigned char*>(device) + words * sizeof(Word);
  queue->parallel_for(sycl::range<1>(std::max(words, rest)), [=](sycl::id<1> index) {
    const auto i = index[0];
    if (i < words) {
      target[i] = source[i];
    }
    if (i < rest) {
      targetRest[i] = sourceRest[i];
    }
  });
}

} // namespace seissol::kernels::device::aux::hostcopy
