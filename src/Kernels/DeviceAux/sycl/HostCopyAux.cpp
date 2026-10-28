// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/HostCopyAux.h"

#include <cstddef>
#include <cstdint>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::hostcopy {

namespace {
void copy(void* target, const void* source, std::size_t bytes, void* stream) {
  if (bytes == 0) {
    return;
  }
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  auto* to = static_cast<unsigned char*>(target);
  const auto* from = static_cast<const unsigned char*>(source);
  queue->parallel_for(sycl::range<1>(bytes),
                      [=](sycl::id<1> index) { to[index[0]] = from[index[0]]; });
}
} // namespace

void copyToHost(void* hostOnDevice, const void* device, std::size_t bytes, void* stream) {
  copy(hostOnDevice, device, bytes, stream);
}

void copyFromHost(void* device, const void* hostOnDevice, std::size_t bytes, void* stream) {
  copy(device, hostOnDevice, bytes, stream);
}

void waitHost(const std::uint32_t* counterOnDevice, std::uint32_t value, void* stream) {
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  queue->single_task([=]() {
    sycl::atomic_ref<std::uint32_t,
                     sycl::memory_order::acquire,
                     sycl::memory_scope::system,
                     sycl::access::address_space::global_space>
        counter(*const_cast<std::uint32_t*>(counterOnDevice));
    while (counter.load() < value) {
    }
  });
}

void signalHost(std::uint32_t* counterOnDevice, std::uint32_t value, void* stream) {
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  queue->single_task([=]() {
    // what the queue has written before reaches the host first
    sycl::atomic_ref<std::uint32_t,
                     sycl::memory_order::release,
                     sycl::memory_scope::system,
                     sycl::access::address_space::global_space>(*counterOnDevice)
        .store(value);
  });
}

} // namespace seissol::kernels::device::aux::hostcopy
