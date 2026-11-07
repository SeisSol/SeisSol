// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/IntegralAux.h"

#include <cstddef>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::integral {

namespace {
template <typename T>
void launch(const T* const* source,
            double* const* target,
            std::size_t size,
            std::size_t count,
            void* stream) {
  if (count == 0 || size == 0) {
    return;
  }
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  queue->parallel_for(sycl::range<2>(count, size), [=](sycl::item<2> item) {
    const auto block = item.get_id(0);
    const auto i = item.get_id(1);
    target[block][i] += static_cast<double>(source[block][i]);
  });
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
