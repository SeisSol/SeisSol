// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "SyncExecutor.h"

#include <Parallel/Host/CpuExecutor.h>
#include <functional>
#include <memory>
namespace seissol::parallel::host {

struct SyncTask : public Task {
  void wait() override {}
};

void SyncExecutor::start(const std::function<void(CpuExecutor&)>& continuation,
                         const Pinning* pinning) {
  std::invoke(continuation, *this);
}
void SyncExecutor::wait() {}
std::shared_ptr<Task> SyncExecutor::add(int priority,
                                        std::size_t size,
                                        const std::function<void(std::size_t)>& function,
                                        const std::vector<std::shared_ptr<Task>>& pollList) {
  for (const auto& task : pollList) {
    task->wait();
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t i = 0; i < size; ++i) {
    std::invoke(function, i);
  }
  return std::make_shared<SyncTask>();
}

} // namespace seissol::parallel::host
