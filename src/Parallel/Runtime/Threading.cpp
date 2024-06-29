#include "Threading.hpp"
#include <atomic>
#include <omp.h>

namespace seissol::parallel::runtime {

void ThreadingRuntime::add(SimpleTask&& task) {
  std::lock_guard tasksLock(tasksMutex);
  tasks.emplace_back(std::move(task));
}

// emulate a de-facto pthread-like runtime; until all OMP features are implemented properly
void ThreadingRuntime::run() {
  std::vector<SimpleTask> nextTasks;
  std::vector<std::atomic<int>> done;
#pragma omp parallel
  {
    while (true) {
#pragma omp barrier
#pragma omp single
      {
        std::lock_guard tasksLock(tasksMutex);
        nextTasks.resize(0);
        done.resize(0);

        for (auto it = tasks.begin(); it != tasks.end(); ++it) {
          if (std::invoke(it->ready)) {
            nextTasks.push_back(*it);
            done.push_back(omp_get_num_threads());
            it = tasks.erase(it);
          }
        }
      }
      for (std::size_t j = 0; j < nextTasks.size(); ++j) {
        auto& task = nextTasks[j];
#pragma omp for nowait
        for (std::size_t i = 0; i < task.size; ++i) {
          std::invoke(task.function, i);
        }
        auto prev = done[j].fetch_sub(1);
        if (prev == 1) {
          std::invoke(task.completion);
        }
      }
    }
  }
}

} // namespace seissol::parallel::runtime
