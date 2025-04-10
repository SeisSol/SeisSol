// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_LOCALBASE_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_LOCALBASE_H_

#include "Common/Constants.h"
#include "generated_code/kernel.h"
#include <Kernels/Local.h>
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

namespace seissol::kernels::solver::linearck {

class Local : public LocalKernel {
  public:
  void setGlobalData(const CompoundGlobalData& compound) override;
  void computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                       LocalData& data,
                       LocalTmp& tmp,
                       const CellMaterialData* materialData,
                       const CellBoundaryMapping (*cellBoundaryMapping)[4],
                       double time,
                       double timeStepWidth) override;

  void computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                              ConditionalMaterialTable& materialTable,
                              ConditionalIndicesTable& indicesTable,
                              kernels::LocalData::Loader& loader,
                              LocalTmp& tmp,
                              double timeStepWidth,
                              seissol::parallel::runtime::StreamRuntime& runtime) override;

  void evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                      ConditionalIndicesTable& indicesTable,
                                      kernels::LocalData::Loader& loader,
                                      seissol::initializer::Layer& layer,
                                      seissol::initializer::LTS& lts,
                                      double time,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsIntegral(const FaceType faceTypes[4],
                     unsigned int& nonZeroFlops,
                     unsigned int& hardwareFlops) override;

  unsigned bytesIntegral() override;

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
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCALBASE_H_
