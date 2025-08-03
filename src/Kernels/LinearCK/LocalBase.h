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
#include "GeneratedCode/kernel.h"
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

template<typename Cfg>
class Local : public LocalKernel<Cfg> {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I<Cfg>::size()],
                       LTS::Ref& data,
                       LocalTmp<Cfg>& tmp,
                       const CellMaterialData* materialData,
                       const CellBoundaryMapping (*cellBoundaryMapping)[4],
                       double time,
                       double timeStepWidth) override;

  void computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                              ConditionalMaterialTable& materialTable,
                              ConditionalIndicesTable& indicesTable,
                              double timeStepWidth,
                              seissol::parallel::runtime::StreamRuntime& runtime) override;

  void evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                      ConditionalIndicesTable& indicesTable,
                                      LTS::Layer& layer,
                                      double time,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsIntegral(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                     std::uint64_t& nonZeroFlops,
                     std::uint64_t& hardwareFlops) override;

  std::uint64_t bytesIntegral() override;

  protected:
  kernel::volume<Cfg> m_volumeKernelPrototype;
  kernel::localFlux<Cfg> m_localFluxKernelPrototype;
  kernel::localFluxNodal<Cfg> m_nodalLfKrnlPrototype;

  kernel::projectToNodalBoundary<Cfg> m_projectKrnlPrototype;
  kernel::projectToNodalBoundaryRotated<Cfg> m_projectRotatedKrnlPrototype;

  kernels::DirichletBoundary<Cfg> dirichletBoundary;

#ifdef ACL_DEVICE
  kernel::gpu_volume<Cfg> deviceVolumeKernelPrototype;
  kernel::gpu_localFlux<Cfg> deviceLocalFluxKernelPrototype;
  kernel::gpu_localFluxNodal<Cfg> deviceNodalLfKrnlPrototype;
  kernel::gpu_projectToNodalBoundaryRotated<Cfg> deviceProjectRotatedKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCALBASE_H_
