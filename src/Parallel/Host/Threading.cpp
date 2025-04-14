// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Threading.h"
#include <Parallel/Host/CpuExecutor.h>
#include <atomic>
#include <omp.h>
#include <unordered_set>

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
  std::vector<std::shared_ptr<Task>> dependencies;
  std::function<bool()> ready;
  std::function<void(std::size_t)> function;
  std::function<void()> completion;
  std::shared_ptr<SimpleEvent> result;
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
       SimpleTask{pollList,
                  [=]() {
                    for (const auto& event : pollList) {
                      if (!dynamic_cast<SimpleEvent*>(event.get())->poll()) {
                        return false;
                      }
                    }
                    return true;
                  },
                  function,
                  [=]() { returnEvent->record(); },
                  returnEvent,
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

struct ExecuteStatus {
  std::vector<SimpleTask> nextTasks;
  std::vector<std::unique_ptr<std::atomic<int>>> taskDone;
  std::atomic<bool> allEmpty{true};

  void fill(std::map<int, std::list<std::shared_ptr<SimpleTask>>>& tasksMap,
            std::mutex& tasksMutex,
            std::unordered_set<void*>& taskset) {
    /*while(this->workDone.load() != 0) {
      // spin wait
    }*/

    std::lock_guard tasksLock(tasksMutex);

    this->nextTasks.resize(0);
    this->taskDone.resize(0);

    this->allEmpty = true;

    std::unordered_set<void*> newtaskset;

    for (auto& [_, tasks] : tasksMap) {
      if (!tasks.empty()) {
        this->allEmpty = false;
      }
      auto counter = 0;
      for (auto it = tasks.begin(); it != tasks.end(); ++it) {
        bool ready = true;
        for (const auto& event : it->get()->dependencies) {
          if (!dynamic_cast<SimpleEvent*>(event.get())->poll() &&
              taskset.find(event.get()) == taskset.end()) {
            ready = false;
            break;
          }
        }
        if (ready) {
          this->nextTasks.emplace_back(**it);
          newtaskset.emplace(it->get()->result.get());
          this->taskDone.emplace_back(std::make_unique<std::atomic<int>>(omp_get_num_threads()));
          it = tasks.erase(it);
        }
        ++counter;
        if (counter >= 30 && !nextTasks.empty()) {
          break;
        }
      }
    }
    taskset.insert(newtaskset.begin(), newtaskset.end());
  }
};

// emulate a de-facto pthread-like runtime; until all OMP features are implemented properly
void ThreadStackExecutor::run() {
  std::vector<ExecuteStatus> executes(16);
  std::unordered_set<void*> tasks;

#pragma omp parallel shared(executes, tasks)
  {
#pragma omp single
    {
      for (auto& execute : executes) {
        execute.fill(tasksMap, tasksMutex, tasks);
      }
    }
    int index = executes.size() - 1;
    while (running.load()) {
      index = (index + 1) % executes.size();
      if (index == 0) {
#pragma omp barrier
#pragma omp single
        {
          for (auto& execute : executes) {
            execute.fill(tasksMap, tasksMutex, tasks);
          }
        }
        if (executes[0].allEmpty.load() && completing.load()) {
          break;
        }
      }
      std::size_t rest = 0;
      for (std::size_t j = 0; j < executes[index].nextTasks.size(); ++j) {
        auto& task = executes[index].nextTasks[j];

        /*const auto perthread = task.size / omp_get_num_threads();
        const auto thisthread = (task.size % omp_get_num_threads() > (omp_get_thread_num() + rest) %
        omp_get_num_threads()) ? perthread + 1 : perthread; const auto prestart = perthread *
        omp_get_thread_num(); const auto truestart = prestart + std::min(task.size %
        omp_get_num_threads(), (omp_get_thread_num() + rest) % omp_get_num_threads()); const auto
        trueend = truestart + thisthread; if (task.size < trueend) { logInfo() <<
        omp_get_thread_num() << truestart << trueend << perthread << thisthread << task.size;
        }
        for (std::size_t i = truestart; i < trueend; ++i) {
          std::invoke(task.function, i);
        }
        //rest = (rest + task.size) % omp_get_num_threads();
        */

        while (!std::invoke(task.ready)) {
        }

#pragma omp for nowait schedule(guided)
        for (std::size_t i = 0; i < task.size; ++i) {
          std::invoke(task.function, i);
        }
        auto prev = executes[index].taskDone[j]->fetch_sub(1);
        if (prev == 1) {
          std::invoke(task.completion);
        }
      }
      // #pragma omp barrier
      // #pragma omp barrier
      /*auto prev = executes[index].workDone.fetch_sub(1);
      if (prev == 1) {
        executes[index].fill(tasksMap, tasksMutex, tasks);
      }*/
    }
  }
  completed.store(true);
}

} // namespace seissol::parallel::host
