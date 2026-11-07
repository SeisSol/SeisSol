// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/ClockAux.h"

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
#endif

namespace seissol::kernels::device::aux::clock {

namespace {
// The host functions after a kernel on the stream run once it is done, and see what it wrote then;
// a fence would only delay it, e.g. behind a large copy to the host.
__global__ void setClockKernel(double* device, double* host, double time) {
  *device = time;
  *host = time;
}

__global__ void advanceClockKernel(double* device, double* host, double timeStepSize) {
  const double time = *device + timeStepSize;
  *device = time;
  *host = time;
}
} // namespace

void set(double* device, double* host, double time, void* stream) {
  setClockKernel<<<1, 1, 0, reinterpret_cast<StreamT>(stream)>>>(device, host, time);
}

void advance(double* device, double* host, double timeStepSize, void* stream) {
  advanceClockKernel<<<1, 1, 0, reinterpret_cast<StreamT>(stream)>>>(device, host, timeStepSize);
}

} // namespace seissol::kernels::device::aux::clock
