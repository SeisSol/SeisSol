// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Touch.h"

#include "GeneratedCode/tensor.h"
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <cstddef>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::kernels {

template<typename RealT>
void touchBuffersDerivatives(RealT** buffers, RealT** derivatives, unsigned numberOfCells) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < numberOfCells; ++cell) {
    // touch buffers
    RealT* buffer = buffers[cell];
    if (buffer != nullptr) {
      for (std::size_t dof = 0; dof < tensor::Q<Cfg>::size(); ++dof) {
        // zero time integration buffers
        buffer[dof] = static_cast<RealT>(0);
      }
    }

    // touch derivatives
    RealT* derivative = derivatives[cell];
    if (derivative != nullptr) {
      for (std::size_t dof = 0; dof < seissol::kernels::Solver<Cfg>::DerivativesSize; ++dof) {
        derivative[dof] = static_cast<RealT>(0);
      }
    }
  }
}

template void touchBuffersDerivatives(double** buffers, double** derivatives, unsigned numberOfCells);
template void touchBuffersDerivatives(float** buffers, float** derivatives, unsigned numberOfCells);

template<typename RealT>
void fillWithStuff(RealT* buffer, unsigned nValues, [[maybe_unused]] bool onDevice) {
  // No RealT point for these numbers. Should be just something != 0 and != NaN and != Inf
  const auto stuff = [](unsigned n) {
    return static_cast<RealT>((214013.0 * n + 2531011.0) / 16777216.0);
  };
#ifdef ACL_DEVICE
  if (onDevice) {
    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();

    device::DeviceInstance::getInstance().algorithms.fillArray<RealT>(
        buffer, static_cast<RealT>(2531011.0 / 65536.0), nValues, stream);

    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
    return;
  }
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned n = 0; n < nValues; ++n) {
    buffer[n] = stuff(n);
  }
}

template void fillWithStuff(double* buffer, unsigned nValues, [[maybe_unused]] bool onDevice);
template void fillWithStuff(float* buffer, unsigned nValues, [[maybe_unused]] bool onDevice);

} // namespace seissol::kernels
