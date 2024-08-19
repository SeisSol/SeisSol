// Copyright (C) 2013-2023 SeisSol group
// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "Touch.h"

#include "generated_code/tensor.h"
#include <utils/logger.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::kernels {

void touchBuffersDerivatives(real** buffers, real** derivatives, unsigned numberOfCells) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    // touch buffers
    // logInfo() << "Reached here 1";
    real* buffer = buffers[cell];
    // logInfo() << "Reached here 2";
    if (buffer != NULL) {
    // logInfo() << "Reached here 3";
      for (unsigned dof = 0; dof < tensor::Q::size(); ++dof) {
        // zero time integration buffers
      // logInfo() << "Reached here 4";
        buffer[dof] = (real)0;
    // logInfo() << "Reached here 5";
      }
    }
    // logInfo() << "Reached here 6";

    // touch derivatives
    real* derivative = derivatives[cell];
    // logInfo() << "Reached here 7";
    if (derivative != NULL) {
    // logInfo() << "Reached here 8";
      for (unsigned dof = 0; dof < yateto::computeFamilySize<tensor::dQ>(); ++dof) {
    // logInfo() << "Reached here 9";
        derivative[dof] = (real)0;
    // logInfo() << "Reached here 10";
      }
    }
  }
}

void fillWithStuff(real* buffer, unsigned nValues, [[maybe_unused]] bool onDevice) {
  // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
  const auto stuff = [](unsigned n) { return static_cast<real>((214013 * n + 2531011) / 65536); };
#ifdef ACL_DEVICE
  if (onDevice) {
    void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();

    device::DeviceInstance::getInstance().algorithms.fillArray<real>(
        buffer, static_cast<real>(2531011.0 / 65536), nValues, stream);

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
