// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "HelperThread.h"

#include "Parallel/Pin.h"

#include <functional>
#include <utility>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel {

HelperThread::HelperThread(std::function<bool()> function, const Pinning* pinning)
    : function_(std::move(function)), pinning_(pinning), isFinished_(false), shouldReset_(false) {}

[[nodiscard]] bool HelperThread::finished() const { return isFinished_.load(); }

void HelperThread::stop() {
  // Send signal to the helper thread to finish and wait
  shouldReset_.store(true);
  if (thread_.joinable()) {
    thread_.join();
  }
}

void HelperThread::start() {
  if (!thread_.joinable()) {
    // Reset flags and reset ghost clusters
    shouldReset_.store(false);
    isFinished_.store(false);

    // Start a new helper thread.
    // Note: Easier than keeping one alive, and not that expensive.
    thread_ = std::thread([this]() {
#ifdef ACL_DEVICE
      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      device.api->setDevice(0);
#endif // ACL_DEVICE
      // Pin this thread to the last core
      // We compute the mask outside the thread because otherwise
      // it confuses profilers and debuggers.
      pinning_->pinToFreeCPUs();
      while (!shouldReset_.load() && !isFinished_.load()) {
        isFinished_.store(std::invoke(function_));
      }
    });
  }
}

void HelperThread::restart() {
  stop();
  start();
}

HelperThread::~HelperThread() {
  if (thread_.joinable()) {
    thread_.join();
  }
}

} // namespace seissol::parallel
