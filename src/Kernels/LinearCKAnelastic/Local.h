// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCAL_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCAL_H_

#include "GeneratedCode/kernel.h"
#include "Kernels/Interface.h"
#include "Kernels/Local.h"
#include "Physics/InitialField.h"

#include <memory>

namespace seissol::kernels::solver::linearckanelastic {
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
  kernel::volumeExt m_volumeKernelPrototype;
  kernel::localFluxExt m_localFluxKernelPrototype;
  kernel::local m_localKernelPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_volumeExt deviceVolumeKernelPrototype;
  kernel::gpu_localFluxExt deviceLocalFluxKernelPrototype;
  kernel::gpu_local deviceLocalKernelPrototype;
  kernel::gpu_fluxLocalAll deviceFluxLocalAllKernelPrototype;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCAL_H_
