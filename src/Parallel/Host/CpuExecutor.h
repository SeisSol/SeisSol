// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_HOST_CPUEXECUTOR_H_
#define SEISSOL_SRC_PARALLEL_HOST_CPUEXECUTOR_H_

#include <Parallel/Pin.h>
#include <atomic>
#include <functional>
#include <memory>
namespace seissol::parallel::host {

class Task {
  public:
  virtual ~Task() = default;
  virtual void wait() = 0;
};

class SimpleEvent : public Task {
  private:
  std::shared_ptr<std::atomic<bool>> data;

  public:
  SimpleEvent() : data(std::make_shared<std::atomic<bool>>(false)) {}
  void init() { data->store(false); }
  void record() { data->store(true); }
  bool poll() { return data->load(); }

  void wait() override {
    while (!poll()) {
      // TODO: yield?
    }
  }
};

class CpuExecutor {
  public:
  virtual ~CpuExecutor() = default;
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
#endif // SEISSOL_SRC_PARALLEL_HOST_CPUEXECUTOR_H_
