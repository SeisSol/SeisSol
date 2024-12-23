// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 * http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Alexander Heinecke (Intel Corp.)
 */

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
