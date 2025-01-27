// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_NEIGHBORBASE_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_NEIGHBORBASE_H_

#include "Common/Constants.h"
#include "generated_code/kernel.h"
#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels {

class NeighborBase {
  protected:
  static void checkGlobalData(const GlobalData* global, size_t alignment);
  kernel::neighboringFlux m_nfKrnlPrototype;
  dynamicRupture::kernel::nodalFlux m_drKrnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_NEIGHBORBASE_H_
