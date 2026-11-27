// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "HelperThread.h"

#include "Core/Parallel/Pin.h"

#include <functional>
#include <utility>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::parallel {

HelperThread::HelperThread(std::function<bool()> function, const Pinning* pinning)
    : function_(std::move(function)), isFinished_(false), shouldReset_(false),
      freeMask_(pinning->getFreeCPUsMask()) {}

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
    thread_ = std::thread([this] {
#ifdef ACL_DEVICE
      // the current device is thread-local; the Device library remembers the one this rank is
      // bound to
      const device::DeviceInstance& device = device::DeviceInstance::instance();
      device.api().setDevice(device.api().getDeviceId());
#endif // ACL_DEVICE
      // Pin this thread to the free CPUs. The mask comes from outside the thread: computing it is
      // a collective on the node, and it would confuse profilers and debuggers.
      Pinning::pinCurrentThread(freeMask_);
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
