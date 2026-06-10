// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Local.h"

#include <memory>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop
#include "Physics/InitialField.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels::solver::linearck {

class Local : public LocalKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeIntegral(real* timeIntegratedDoFs,
                       LTS::Ref& data,
                       LocalTmp& tmp,
                       double time,
                       double timeStepWidth) override;

  void computeBatchedIntegral(recording::ConditionalPointersToRealsTable& dataTable,
                              recording::ConditionalMaterialTable& materialTable,
                              recording::ConditionalIndicesTable& indicesTable,
                              double timeStepWidth,
                              seissol::parallel::runtime::StreamRuntime& runtime) override;

  void evaluateBatchedTimeDependentBc(recording::ConditionalPointersToRealsTable& dataTable,
                                      recording::ConditionalIndicesTable& indicesTable,
                                      LTS::Layer& layer,
                                      double time,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsIntegral(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                     std::uint64_t& nonZeroFlops,
                     std::uint64_t& hardwareFlops) override;

  std::uint64_t bytesIntegral() override;

  protected:
  kernel::volume volumeKernelPrototype_;
  kernel::localFlux localFluxKernelPrototype_;
  kernel::localFluxNodal nodalLfKrnlPrototype_;

  kernel::projectToNodalBoundary projectKrnlPrototype_;
  kernel::projectToNodalBoundaryRotated projectRotatedKrnlPrototype_;

  kernels::DirichletBoundary dirichletBoundary_;

  kernel::bcFreeSurfaceGravity bcFreeSurfaceGravity_;
  kernel::bcDirichlet bcDirichlet_;

#ifdef ACL_DEVICE
  kernel::gpu_volume deviceVolumeKernelPrototype_;
  kernel::gpu_localFlux deviceLocalFluxKernelPrototype_;
  kernel::gpu_localFluxAll deviceLocalFluxAllKernelPrototype_;
  kernel::gpu_localFluxNodal deviceNodalLfKrnlPrototype_;
  kernel::gpu_projectToNodalBoundaryRotated deviceProjectRotatedKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();

  kernel::gpu_bcFreeSurfaceGravity deviceBCFreeSurfaceGravity_;
  kernel::gpu_bcDirichlet deviceBCDirichlet_;
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_
