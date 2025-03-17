// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_HOST_THREADING_H_
#define SEISSOL_SRC_PARALLEL_HOST_THREADING_H_

#include <Parallel/HelperThread.h>
#include <Parallel/Host/CpuExecutor.h>
#include <atomic>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <mutex>
namespace seissol::parallel::host {

class SimpleTask;

class ThreadStackExecutor : public CpuExecutor {
  public:
  ~ThreadStackExecutor() override = default;
  void start(const std::function<void(CpuExecutor&)>& continuation,
             const Pinning* pinning) override;
  std::shared_ptr<Task> add(int priority,
                            std::size_t size,
                            const std::function<void(std::size_t)>& function,
                            const std::vector<std::shared_ptr<Task>>& pollList) override;
  // std::shared_ptr<Task> fromDevice(void* stream) override;

  private:
  void task(int priority, SimpleTask&& task);
  void run();
  void complete();
  std::atomic<bool> running;
  std::atomic<bool> completed;
  std::mutex tasksMutex;
  std::map<int, std::list<SimpleTask>> tasksMap;
};

} // namespace seissol::parallel::host

#endif // SEISSOL_SRC_PARALLEL_HOST_THREADING_H_
