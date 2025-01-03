#pragma once

#include <Parallel/Pin.h>
#include <atomic>
#include <functional>
#include <thread>
namespace seissol::parallel {

class HelperThread {
  public:
  HelperThread(std::function<bool()> function, const Pinning* pinning);
  void stop();
  void start();
  void restart();
  [[nodiscard]] bool finished() const;

  ~HelperThread();

  private:
  std::function<bool()> function;
  std::thread thread;
  std::atomic<bool> shouldReset;
  std::atomic<bool> isFinished;
  const Pinning* pinning;
};

} // namespace seissol::parallel
