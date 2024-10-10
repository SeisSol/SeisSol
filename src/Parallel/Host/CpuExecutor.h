#pragma once

#include <Parallel/Pin.h>
#include <functional>
#include <memory>
namespace seissol::parallel::host {

class Task {};

class CpuExecutor {
  public:
  virtual void start(const std::function<void(CpuExecutor&)>& continuation,
                     const Pinning* pinning) = 0;
  virtual void wait() = 0;
  virtual std::shared_ptr<Task> add(int priority,
                                    std::size_t size,
                                    const std::function<void(std::size_t)>& function,
                                    const std::vector<std::shared_ptr<Task>>& pollList) = 0;
  // virtual std::shared_ptr<Task> fromDevice(void* stream);
};

} // namespace seissol::parallel::host
