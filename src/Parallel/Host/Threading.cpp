#include "Threading.h"
#include <atomic>
#include <omp.h>

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
};
#endif

class SimpleEvent : public Task {
  private:
  std::shared_ptr<std::atomic<bool>> data;

  public:
  SimpleEvent() : data(std::make_shared<std::atomic<bool>>(false)) {}
  void init() { data->store(false); }
  void record() { data->store(true); }
  bool poll() { return data->load(); }
};

struct SimpleTask {
  std::function<bool()> ready;
  std::function<void(std::size_t)> function;
  std::function<void()> completion;
  std::size_t size;
};

void ThreadStackExecutor::start(const std::function<void(CpuExecutor&)>& continuation,
                                const Pinning* pinning) {
  running.store(true);
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
  tasksMap[priority].emplace_back(std::move(task));
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
                      if (!static_cast<SimpleEvent*>(event.get())->poll()) {
                        return false;
                      }
                    }
                    return true;
                  },
                  function,
                  [=]() { returnEvent->record(); }});
  return returnEvent;
}

/*
std::shared_ptr<Task> ThreadStackExecutor::fromDevice(void* stream) {
  return std::make_shared<FromDeviceEvent>();
}
*/

void ThreadStackExecutor::complete() { completed.store(true); }

// emulate a de-facto pthread-like runtime; until all OMP features are implemented properly
void ThreadStackExecutor::run() {
  std::vector<SimpleTask> nextTasks;
  std::vector<std::atomic<int>> done;
#pragma omp parallel
  {
    while (running.load()) {
      bool allEmpty = true;
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
            if (std::invoke(it->ready)) {
              nextTasks.push_back(*it);
              done.push_back(omp_get_num_threads());
              it = tasks.erase(it);
            }
          }
        }
      }
      if (allEmpty && completed.load()) {
        break;
      }
      for (std::size_t j = 0; j < nextTasks.size(); ++j) {
        auto& task = nextTasks[j];
#pragma omp for nowait schedule(dynamic, 1)
        for (std::size_t i = 0; i < task.size; ++i) {
          std::invoke(task.function, i);
        }
        auto prev = done[j].fetch_sub(1);
        if (prev == 1) {
          std::invoke(task.completion);
        }
      }
#pragma omp barrier
    }
  }
}

} // namespace seissol::parallel::host
