// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ClusterClock.h"

#include "Parallel/Runtime/Stream.h"

#ifdef ACL_DEVICE
#include "Kernels/DeviceAux/ClockAux.h"

#include <Device/AbstractAPI.h>
#include <Device/device.h>
#endif

namespace seissol::solver {

#ifdef ACL_DEVICE

namespace {
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::instance(); }
} // namespace

ClusterClock::ClusterClock()
    : host_(static_cast<double*>(deviceInstance().api().allocPinnedMem(
          sizeof(double), false, device::Destination::CurrentDevice))),
      device_(static_cast<double*>(deviceInstance().api().allocGlobMem(sizeof(double)))),
      hostOnDevice_(static_cast<double*>(deviceInstance().api().devicePointer(host_))) {
  *host_ = 0;
  const double zero = 0;
  deviceInstance().api().copyTo(device_, &zero, sizeof(double));
}

ClusterClock::~ClusterClock() { dispose(); }

void ClusterClock::dispose() {
  if (host_ != nullptr) {
    deviceInstance().api().freeGlobMem(device_);
    deviceInstance().api().freePinnedMem(host_);
    host_ = nullptr;
    device_ = nullptr;
    hostOnDevice_ = nullptr;
  }
}

void ClusterClock::set(double time, parallel::runtime::StreamRuntime& runtime) {
  kernels::device::aux::clock::set(device_, hostOnDevice_, time, runtime.stream());
}

void ClusterClock::advance(double timeStepSize, parallel::runtime::StreamRuntime& runtime) {
  kernels::device::aux::clock::advance(device_, hostOnDevice_, timeStepSize, runtime.stream());
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
