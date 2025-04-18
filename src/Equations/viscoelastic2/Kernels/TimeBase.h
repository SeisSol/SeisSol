// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_TIMEBASE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_TIMEBASE_H_

#include "generated_code/kernel.h"

namespace seissol::kernels {
class TimeBase {
  protected:
  kernel::derivative m_krnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_TIMEBASE_H_
