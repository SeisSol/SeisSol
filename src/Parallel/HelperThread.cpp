// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "HelperThread.h"

#include <utility>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::parallel {

HelperThread::HelperThread(std::function<bool()> function, const Pinning* pinning)
    : function(std::move(function)), pinning(pinning), isFinished(false), shouldReset(false) {}

[[nodiscard]] bool HelperThread::finished() const { return isFinished.load(); }

void HelperThread::stop() {
  // Send signal to the helper thread to finish and wait
  shouldReset.store(true);
  if (thread.joinable()) {
    thread.join();
  }
}

void HelperThread::start() {
  if (!thread.joinable()) {
    // Reset flags and reset ghost clusters
    shouldReset.store(false);
    isFinished.store(false);

    // Start a new helper thread.
    // Note: Easier than keeping one alive, and not that expensive.
    thread = std::thread([this]() {
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      device.api->setDevice(0);
#endif // ACL_DEVICE
      // Pin this thread to the last core
      // We compute the mask outside the thread because otherwise
      // it confuses profilers and debuggers.
      pinning->pinToFreeCPUs();
      while (!shouldReset.load() && !isFinished.load()) {
        isFinished.store(std::invoke(function));
      }
    });
  }
}

void HelperThread::restart() {
  stop();
  start();
}

HelperThread::~HelperThread() {
  if (thread.joinable()) {
    thread.join();
  }
}

} // namespace seissol::parallel
