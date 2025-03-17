// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_HELPERTHREAD_H_
#define SEISSOL_SRC_PARALLEL_HELPERTHREAD_H_

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
#endif // SEISSOL_SRC_PARALLEL_HELPERTHREAD_H_
