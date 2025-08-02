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

void touchBuffersDerivatives(real** buffers, real** derivatives, unsigned numberOfCells) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < numberOfCells; ++cell) {
    // touch buffers
    real* buffer = buffers[cell];
    if (buffer != nullptr) {
      for (std::size_t dof = 0; dof < tensor::Q<Cfg>::size(); ++dof) {
        // zero time integration buffers
        buffer[dof] = static_cast<real>(0);
      }
    }

    // touch derivatives
    real* derivative = derivatives[cell];
    if (derivative != nullptr) {
      for (std::size_t dof = 0; dof < seissol::kernels::Solver::DerivativesSize; ++dof) {
        derivative[dof] = static_cast<real>(0);
      }
    }
  }
}

void fillWithStuff(real* buffer, unsigned nValues, [[maybe_unused]] bool onDevice) {
  // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
  const auto stuff = [](unsigned n) {
    return static_cast<real>((214013.0 * n + 2531011.0) / 16777216.0);
  };
#ifdef ACL_DEVICE
  if (onDevice) {
    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();

    device::DeviceInstance::getInstance().algorithms.fillArray<real>(
        buffer, static_cast<real>(2531011.0 / 65536.0), nValues, stream);

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

} // namespace seissol::kernels
