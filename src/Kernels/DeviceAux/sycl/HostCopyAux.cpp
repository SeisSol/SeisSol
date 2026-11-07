// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/HostCopyAux.h"

#include <cstddef>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::hostcopy {

void copyToHost(void* hostOnDevice, const void* device, std::size_t bytes, void* stream) {
  if (bytes == 0) {
    return;
  }
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  auto* target = static_cast<unsigned char*>(hostOnDevice);
  const auto* source = static_cast<const unsigned char*>(device);
  queue->parallel_for(sycl::range<1>(bytes),
                      [=](sycl::id<1> index) { target[index[0]] = source[index[0]]; });
}

} // namespace seissol::kernels::device::aux::hostcopy
