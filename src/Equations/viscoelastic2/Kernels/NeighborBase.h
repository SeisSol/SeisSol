// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_

#include "generated_code/kernel.h"

namespace seissol::kernels {
class NeighborBase {
  protected:
  kernel::neighborFluxExt m_nfKrnlPrototype;
  kernel::neighbor m_nKrnlPrototype;
  dynamicRupture::kernel::nodalFlux m_drKrnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt deviceNfKrnlPrototype;
  kernel::gpu_neighbor deviceNKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_
