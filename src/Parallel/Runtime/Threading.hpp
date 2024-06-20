#ifndef SEISSOL_PARALLEL_RUNTIME_THREADING_HPP
#define SEISSOL_PARALLEL_RUNTIME_THREADING_HPP

#include <functional>
namespace seissol::parallel::runtime {

class ThreadingRuntime {
  public:
  void run();
  void add();
  void addFor();

  private:
  std::vector<std::function<void(std::size_t)>> functions;
};

} // namespace seissol::parallel::runtime

#endif
