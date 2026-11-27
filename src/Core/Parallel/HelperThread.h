// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_CORE_PARALLEL_HELPERTHREAD_H_
#define SEISSOL_SRC_CORE_PARALLEL_HELPERTHREAD_H_

#include "Core/Parallel/Pin.h"

#include <async/as/Pin.h>
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

  HelperThread(const HelperThread&) = delete;
  auto operator=(const HelperThread&) -> HelperThread& = delete;

  HelperThread(HelperThread&&) = delete;
  auto operator=(HelperThread&&) -> HelperThread& = delete;

  private:
  std::function<bool()> function_;
  std::atomic<bool> isFinished_;
  std::atomic<bool> shouldReset_;
  std::thread thread_;
  // the free CPUs, computed by the thread that creates the helper: the computation is a collective
  // on the processes of the node, which the helper thread must not call (it starts anew in each
  // epoch)
  async::as::CpuMask freeMask_;
};

} // namespace seissol::parallel
#endif // SEISSOL_SRC_CORE_PARALLEL_HELPERTHREAD_H_
