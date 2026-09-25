// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ClusterClock.h"

#include "Parallel/Runtime/Stream.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

#ifdef ACL_DEVICE

namespace {
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::getInstance(); }
} // namespace

ClusterClock::ClusterClock() {
  host_ = static_cast<double*>(deviceInstance().api->allocPinnedMem(sizeof(double)));
  device_ = static_cast<double*>(deviceInstance().api->allocGlobMem(sizeof(double)));
  deviceStep_ = static_cast<double*>(deviceInstance().api->allocGlobMem(sizeof(double)));
  stepTable_ = static_cast<const double**>(deviceInstance().api->allocGlobMem(sizeof(double*)));
  clockTable_ = static_cast<double**>(deviceInstance().api->allocGlobMem(sizeof(double*)));

  *host_ = 0;
  const double zero = 0;
  deviceInstance().api->copyTo(device_, &zero, sizeof(double));
  deviceInstance().api->copyTo(deviceStep_, &zero, sizeof(double));
  deviceInstance().api->copyTo(stepTable_, &deviceStep_, sizeof(double*));
  deviceInstance().api->copyTo(clockTable_, &device_, sizeof(double*));
}

ClusterClock::~ClusterClock() { dispose(); }

void ClusterClock::dispose() {
  if (host_ != nullptr) {
    deviceInstance().api->freeGlobMem(clockTable_);
    deviceInstance().api->freeGlobMem(stepTable_);
    deviceInstance().api->freeGlobMem(deviceStep_);
    deviceInstance().api->freeGlobMem(device_);
    deviceInstance().api->freePinnedMem(host_);
    host_ = nullptr;
    device_ = nullptr;
  }
}

void ClusterClock::set(double time, parallel::runtime::StreamRuntime& runtime) {
  deviceInstance().algorithms.fillArray(device_, time, 1, runtime.stream());
  deviceInstance().api->copyFromAsync(host_, device_, sizeof(double), runtime.stream());
}

void ClusterClock::advance(double timeStepSize, parallel::runtime::StreamRuntime& runtime) {
  deviceInstance().algorithms.fillArray(deviceStep_, timeStepSize, 1, runtime.stream());
  deviceInstance().algorithms.accumulateBatchedData(
      stepTable_, clockTable_, 1, 1, runtime.stream());
  deviceInstance().api->copyFromAsync(host_, device_, sizeof(double), runtime.stream());
}

#else

ClusterClock::ClusterClock() : host_(new double(0)) {}

ClusterClock::~ClusterClock() { dispose(); }

void ClusterClock::dispose() {
  delete host_;
  host_ = nullptr;
}

void ClusterClock::set(double time, parallel::runtime::StreamRuntime& /*runtime*/) {
  *host_ = time;
}

void ClusterClock::advance(double timeStepSize, parallel::runtime::StreamRuntime& /*runtime*/) {
  *host_ += timeStepSize;
}

#endif

} // namespace seissol::solver
