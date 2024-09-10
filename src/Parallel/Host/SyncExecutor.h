#pragma once

#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Pin.h>
#include <functional>
#include <memory>
namespace seissol::parallel::host {

class SyncExecutor : public CpuExecutor {
  public:
  void start(const std::function<void(CpuExecutor&)>& continuation,
             const Pinning* pinning) override;
  void wait() override;
  std::shared_ptr<Task> add(int priority,
                            std::size_t size,
                            const std::function<void(std::size_t)>& function,
                            const std::vector<std::shared_ptr<Task>>& pollList) override;
  // std::shared_ptr<Task> fromDevice(void* stream) override;
};

} // namespace seissol::parallel::host
