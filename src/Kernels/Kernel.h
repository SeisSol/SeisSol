// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_KERNEL_H_
#define SEISSOL_SRC_KERNELS_KERNEL_H_

#include <Common/Executor.h>
#include <Initializer/Typedefs.h>
#include <Parallel/Runtime/Stream.h>

namespace seissol::kernels {

class Kernel {
  public:
  virtual void setGlobalData(const CompoundGlobalData& global) {}
  virtual ~Kernel() = default;

#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels
#endif // SEISSOL_SRC_KERNELS_KERNEL_H_
