#ifndef SEISSOL_PARALLEL_RUNTIME_THREADING_HPP
#define SEISSOL_PARALLEL_RUNTIME_THREADING_HPP

#include <Parallel/HelperThread.hpp>
#include <atomic>
#include <functional>
#include <list>
#include <map>
#include <mutex>
namespace seissol::parallel::runtime {

struct SimpleTask {
  std::function<bool()> ready;
  std::function<void(std::size_t)> function;
  std::function<void()> completion;
  std::size_t size;
};

class ThreadingRuntime {
  public:
  template <typename F>
  void start(F&& continueExecution) {
    running.store(true);
    HelperThread thread([this, continueExecution] {
      std::invoke(continueExecution);
      running.store(false);
    });
    thread.start();
    run();
  }

  void abort();
  void add(int priority, SimpleTask&& task);

  private:
  void run();
  std::atomic<bool> running;
  std::mutex tasksMutex;
  std::map<int, std::list<SimpleTask>> tasksMap;
};

} // namespace seissol::parallel::runtime

#endif
