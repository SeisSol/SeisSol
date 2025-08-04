// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCALBASE_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCALBASE_H_

#include "GeneratedCode/kernel.h"
#include "Physics/InitialField.h"
#include <Kernels/Interface.h>
#include <Kernels/Local.h>
#include <memory>

namespace seissol::kernels::solver::linearckanelastic {

template<typename Cfg>
class Local : public LocalKernel<Cfg> {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I<Cfg>::size()],
                       LTS::Ref<Cfg>& data,
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
  kernel::volumeExt<Cfg> m_volumeKernelPrototype;
  kernel::localFluxExt<Cfg> m_localFluxKernelPrototype;
  kernel::local<Cfg> m_localKernelPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_volumeExt<Cfg> deviceVolumeKernelPrototype;
  kernel::gpu_localFluxExt<Cfg> deviceLocalFluxKernelPrototype;
  kernel::gpu_local<Cfg> deviceLocalKernelPrototype;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCALBASE_H_
