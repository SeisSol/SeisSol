// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "OMPTasking.h"

#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Pin.h>
#include <functional>
#include <memory>
#include <omp.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::parallel::host {

class OMPTask : public Task {
  omp_depend_t depobj;

  void wait() override {
#pragma omp taskwait // TODO
  }
};

void OMPTaskingExecutor::start(const std::function<void(CpuExecutor&)>& continuation,
                               const Pinning* pinning) {
#pragma omp parallel
  {
#ifdef ACL_DEVICE
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.api->setDevice(0);
#endif
#pragma omp single
    {
      std::invoke(continuation, *this);
    }
  }
#pragma omp taskwait
}
std::shared_ptr<Task> OMPTaskingExecutor::add(int priority,
                                              std::size_t size,
                                              const std::function<void(std::size_t)>& function,
                                              const std::vector<std::shared_ptr<Task>>& pollList) {
  omp_depend_t depobj;
#pragma omp taskloop priority(priority)
  for (std::size_t i = 0; i < size; ++i) {
    std::invoke(function, i);
  }
  return std::make_shared<OMPTask>();
}

} // namespace seissol::parallel::host
