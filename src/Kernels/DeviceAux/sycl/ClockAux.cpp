// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/ClockAux.h"

#include "Kernels/DeviceAux/TimePointsAux.h"

#include <cstddef>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::clock {

void set(double* device, double* host, double time, void* stream) {
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  queue->single_task([=] {
    *device = time;
    *host = time;
  });
}

void advance(double* device, double* host, double timeStepSize, void* stream) {
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  queue->single_task([=] {
    const double time = *device + timeStepSize;
    *device = time;
    *host = time;
  });
}

} // namespace seissol::kernels::device::aux::clock

namespace seissol::kernels::timepoints {

void timePoints(const double* clock, const TimeOffsets& offsets, double* times, void* stream) {
  if (offsets.count > 0) {
    auto* queue = reinterpret_cast<sycl::queue*>(stream);
    queue->single_task([=] {
      for (std::size_t i = 0; i < offsets.count; ++i) {
        times[i] = *clock + offsets.values[i];
      }
    });
  }
}

} // namespace seissol::kernels::timepoints
