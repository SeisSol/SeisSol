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
  void computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                       LTS::Ref& data,
                       LocalTmp& tmp,
                       const CellMaterialData* materialData,
                       const CellBoundaryMapping (*cellBoundaryMapping)[4],
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

  kernel::gpu_bcFreeSurfaceGravity deviceBCFSG;
  kernel::gpu_bcDirichlet deviceBCDirichlet;
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_
