#include "SyncExecutor.h"

#include <Parallel/Host/CpuExecutor.h>
#include <functional>
#include <memory>
namespace seissol::parallel::host {

void SyncExecutor::start(const std::function<void(CpuExecutor&)>& continuation,
                         const Pinning* pinning) {
  std::invoke(continuation, *this);
}
void SyncExecutor::wait() {}
std::shared_ptr<Task> SyncExecutor::add(int priority,
                                        std::size_t size,
                                        const std::function<void(std::size_t)>& function,
                                        const std::vector<std::shared_ptr<Task>>& pollList) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t i = 0; i < size; ++i) {
    std::invoke(function, i);
  }
  return std::make_shared<Task>();
}

} // namespace seissol::parallel::host
