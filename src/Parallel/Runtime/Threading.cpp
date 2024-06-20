#include "Threading.hpp"

namespace seissol::parallel::runtime {

// emulate a de-facto pthread-like runtime; until all OMP features are implemented properly
void ThreadingRuntime::run() {
#pragma omp parallel
  {
    while (true) {
#pragma omp barrier
      std::vector<std::function<void(std::size_t)>> nextTasks;
      for (const auto& task : nextTasks) {
#pragma omp for nowait
        for (std::size_t i = 0; i < TODO; ++i) {
          task(i);
        }
        // TODO: notify task done
      }
    }
  }
}

} // namespace seissol::parallel::runtime
