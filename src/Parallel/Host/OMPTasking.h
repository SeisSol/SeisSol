#pragma once

#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Pin.h>
#include <functional>
#include <memory>
#include <omp.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::parallel::host {

class OMPTaskingExecutor : public CpuExecutor {
  public:
  void start(const std::function<void(CpuExecutor&)>& continuation,
             const Pinning* pinning) override;
  std::shared_ptr<Task> add(int priority,
                            std::size_t size,
                            const std::function<void(std::size_t)>& function,
                            const std::vector<std::shared_ptr<Task>>& pollList) override;
  // std::shared_ptr<Task> fromDevice(void* stream) override;
};

} // namespace seissol::parallel::host
