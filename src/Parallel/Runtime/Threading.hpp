#ifndef SEISSOL_PARALLEL_RUNTIME_THREADING_HPP
#define SEISSOL_PARALLEL_RUNTIME_THREADING_HPP

#include <functional>
#include <list>
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
  void run();
  void add(SimpleTask&& task);

  private:
  std::mutex tasksMutex;
  std::list<SimpleTask> tasks;
};

} // namespace seissol::parallel::runtime

#endif
