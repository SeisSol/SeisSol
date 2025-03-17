// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Threading.h"
#include <atomic>
#include <omp.h>

#include "Parallel/HostEvent.h"

#include "utils/logger.h"

namespace seissol::parallel::host {

#ifdef ACL_DEVICE
class FromDeviceEvent : public Task {
  private:
  void* event;

  public:
  FromDeviceEvent(void* event) : event(event) {}
  void init() {}
  void record() {}
  bool poll() { return true; }

  void wait() override {
    while (!poll())
      ;
  }
};

class ToDeviceEvent : public Task {
  private:
  seissol::parallel::HostEvent event;

  public:
  ToDeviceEvent(void* stream) { event.streamWait(stream); }

  void init() {}
  void record() { event.complete(); }
  bool poll() { return event.completed(); }

  void wait() override { event.hostWait(); }
};
#endif

struct SimpleTask {
  std::function<bool()> ready;
  std::function<void(std::size_t)> function;
  std::function<void()> completion;
  std::size_t size;
};

void ThreadStackExecutor::start(const std::function<void(CpuExecutor&)>& continuation,
                                const Pinning* pinning) {
  running.store(true);
  completing.store(false);
  completed.store(false);
  auto thread = HelperThread(
      [this, continuation] {
        std::invoke(continuation, *this);
        complete();
        return true;
      },
      pinning);
  thread.start();
  run();
}

void ThreadStackExecutor::task(int priority, SimpleTask&& task) {
  std::lock_guard tasksLock(tasksMutex);
  tasksMap[priority].emplace_back(std::make_shared<SimpleTask>(std::move(task)));
}

std::shared_ptr<Task> ThreadStackExecutor::add(int priority,
                                               std::size_t size,
                                               const std::function<void(std::size_t)>& function,
                                               const std::vector<std::shared_ptr<Task>>& pollList) {
  auto returnEvent = std::make_shared<SimpleEvent>();
  returnEvent->init();
  task(priority,
       SimpleTask{[=]() {
                    for (const auto& event : pollList) {
                      if (!dynamic_cast<SimpleEvent*>(event.get())->poll()) {
                        return false;
                      }
                    }
                    return true;
                  },
                  function,
                  [=]() { returnEvent->record(); },
                  size});
  return returnEvent;
}

/*
std::shared_ptr<Task> ThreadStackExecutor::fromDevice(void* stream) {
  return std::make_shared<FromDeviceEvent>();
}
*/

void ThreadStackExecutor::wait() {
  while (!completed.load()) {
  }
}

void ThreadStackExecutor::complete() { completing.store(true); }

// emulate a de-facto pthread-like runtime; until all OMP features are implemented properly
void ThreadStackExecutor::run() {
  std::vector<SimpleTask> nextTasks;
  std::vector<std::unique_ptr<std::atomic<int>>> done;

  std::atomic<bool> allEmpty{false};
#pragma omp parallel shared(nextTasks, done, allEmpty)
  {
    while (running.load()) {
#pragma omp single
      {
        std::lock_guard tasksLock(tasksMutex);
        nextTasks.resize(0);
        done.resize(0);

        for (auto& [_, tasks] : tasksMap) {
          if (!tasks.empty()) {
            allEmpty = false;
          }
          for (auto it = tasks.begin(); it != tasks.end(); ++it) {
            if (std::invoke(it->get()->ready)) {
              nextTasks.emplace_back(**it);
              done.emplace_back(std::make_unique<std::atomic<int>>(omp_get_num_threads()));
              it = tasks.erase(it);
            }
          }
        }
        allEmpty = done.empty();
      }

      // implicit barrier here
      if (allEmpty && completing.load()) {
        break;
      }
      for (std::size_t j = 0; j < nextTasks.size(); ++j) {
        auto& task = nextTasks[j];
#pragma omp for nowait schedule(dynamic, 4)
        for (std::size_t i = 0; i < task.size; ++i) {
          std::invoke(task.function, i);
        }
        auto prev = done[j]->fetch_sub(1);
        if (prev == 1) {
          std::invoke(task.completion);
        }
      }
#pragma omp barrier
    }
  }
  completed.store(true);
}

} // namespace seissol::parallel::host
