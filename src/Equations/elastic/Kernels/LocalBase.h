// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_LOCALBASE_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_LOCALBASE_H_

#include "Common/Constants.h"
#include "generated_code/kernel.h"
#include <memory>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop
#include "Physics/InitialField.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels {

class LocalBase {
  protected:
  double gravitationalAcceleration;
  static void checkGlobalData(const GlobalData* global, size_t alignment);
  kernel::volume m_volumeKernelPrototype;
  kernel::localFlux m_localFluxKernelPrototype;
  kernel::localFluxNodal m_nodalLfKrnlPrototype;

  kernel::projectToNodalBoundary m_projectKrnlPrototype;
  kernel::projectToNodalBoundaryRotated m_projectRotatedKrnlPrototype;

  kernels::DirichletBoundary dirichletBoundary;

#ifdef ACL_DEVICE
  kernel::gpu_volume deviceVolumeKernelPrototype;
  kernel::gpu_localFlux deviceLocalFluxKernelPrototype;
  kernel::gpu_localFluxNodal deviceNodalLfKrnlPrototype;
  kernel::gpu_projectToNodalBoundaryRotated deviceProjectRotatedKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  const std::vector<std::unique_ptr<physics::InitialField>>* initConds;

  public:
  virtual void setInitConds(decltype(initConds) initConds) { this->initConds = initConds; }

  void setGravitationalAcceleration(double g) { gravitationalAcceleration = g; }

  physics::InitialField* getInitCond(size_t index) {
    const auto& condition = this->initConds->at(index);
    return condition.get();
  }
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_LOCALBASE_H_
