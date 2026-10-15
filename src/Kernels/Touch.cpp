// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Touch.h"

#include "Common/Real.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/SolverSelector.h"

#include <cstddef>
#include <yateto.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels {

template <typename Cfg>
void touchBuffersDerivatives(Real<Cfg>** buffers, Real<Cfg>** derivatives, unsigned numberOfCells) {
  using real = Real<Cfg>;

#pragma omp parallel for schedule(static)
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
      for (std::size_t dof = 0; dof < seissol::kernels::SolverOf<Cfg>::DerivativesSize; ++dof) {
        derivative[dof] = static_cast<real>(0);
      }
    }
  }
}

template <typename RealT>
void fillWithStuff(RealT* buffer, unsigned nValues, [[maybe_unused]] bool onDevice) {
  // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
  const auto stuff = [](unsigned n) {
    return static_cast<RealT>((214013.0 * n + 2531011.0) / 16777216.0);
  };
#ifdef ACL_DEVICE
  if (onDevice) {
    void* stream = device::DeviceInstance::instance().api().getDefaultStream();

    device::DeviceInstance::instance().algorithms().fillArray<RealT>(
        buffer, static_cast<RealT>(2531011.0 / 65536.0), nValues, stream);

    device::DeviceInstance::instance().api().syncDefaultStreamWithHost();
    return;
  }
#endif

#pragma omp parallel for schedule(static)
  for (unsigned n = 0; n < nValues; ++n) {
    buffer[n] = stuff(n);
  }
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template void touchBuffersDerivatives<Cfg>(Real<Cfg>**, Real<Cfg>**, unsigned);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

template void fillWithStuff(float* buffer, unsigned nValues, bool onDevice);
template void fillWithStuff(double* buffer, unsigned nValues, bool onDevice);

} // namespace seissol::kernels
