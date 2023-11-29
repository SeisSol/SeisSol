// Copyright (C) 2013-2023 SeisSol group
// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "Touch.h"

#include <generated_code/tensor.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Parallel/AcceleratorDevice.h"
#endif

namespace seissol::kernels {

void touchBuffersDerivatives(real** buffers, real** derivatives, unsigned numberOfCells) {
#ifdef ACL_DEVICE
  constexpr auto qSize = tensor::Q::size();
  constexpr auto dQSize = yateto::computeFamilySize<tensor::dQ>();
  auto queue = seissol::AcceleratorDevice::getInstance().getSyclDefaultQueue();
  queue
      .parallel_for(sycl::range<1>{numberOfCells},
                    [=](sycl::id<1> idx) {
                      if (real* buffer = buffers[idx[0]]; buffer != NULL) {
                        for (unsigned dof = 0; dof < qSize; ++dof) {
                          buffer[dof] = static_cast<real>(0.0);
                        }
                      }

                      if (real* derivative = derivatives[idx[0]]; derivative != NULL) {
                        for (unsigned dof = 0; dof < dQSize; ++dof) {
                          derivative[dof] = static_cast<real>(0.0);
                        }
                      }
                    })
      .wait();
#else
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    // touch buffers
    real* buffer = buffers[cell];
    if (buffer != NULL) {
      for (unsigned dof = 0; dof < tensor::Q::size(); ++dof) {
        // zero time integration buffers
        buffer[dof] = (real)0;
      }
    }

    // touch derivatives
    real* derivative = derivatives[cell];
    if (derivative != NULL) {
      for (unsigned dof = 0; dof < yateto::computeFamilySize<tensor::dQ>(); ++dof) {
        derivative[dof] = (real)0;
      }
    }
  }
#endif
}

void fillWithStuff(real* buffer, unsigned nValues, [[maybe_unused]] bool onDevice) {
  // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
  auto const stuff = [](unsigned n) { return static_cast<real>((214013 * n + 2531011) / 65536); };
#ifdef ACL_DEVICE
  if (onDevice) {
    auto queue = seissol::AcceleratorDevice::getInstance().getSyclDefaultQueue();
    queue
        .parallel_for(sycl::range<1>{nValues},
                      [=](sycl::id<1> idx) { buffer[idx[0]] = stuff(idx[0]); })
        .wait();
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
