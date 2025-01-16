// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_TIMEBASE_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_TIMEBASE_H_

#include "Common/Constants.h"
#include "generated_code/kernel.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif // ACL_DEVICE

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels {

class TimeBase {
  protected:
  static void checkGlobalData(const GlobalData* global, size_t alignment);
#ifdef USE_STP
  kernel::spaceTimePredictor m_krnlPrototype;
#else
  kernel::derivative m_krnlPrototype;
#endif
  kernel::projectDerivativeToNodalBoundaryRotated projectDerivativeToNodalBoundaryRotated;

  unsigned int m_derivativesOffsets[ConvergenceOrder];

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceKrnlPrototype;
  kernel::gpu_projectDerivativeToNodalBoundaryRotated deviceDerivativeToNodalBoundaryRotated;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  public:
  /**
   * Constructor, which initializes the time kernel.
   **/
  TimeBase();
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_TIMEBASE_H_
