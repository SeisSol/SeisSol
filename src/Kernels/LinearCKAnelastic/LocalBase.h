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
class Local : public LocalKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

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
  kernel::volumeExt m_volumeKernelPrototype;
  kernel::localFluxExt m_localFluxKernelPrototype;
  kernel::local m_localKernelPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_volumeExt deviceVolumeKernelPrototype;
  kernel::gpu_localFluxExt deviceLocalFluxKernelPrototype;
  kernel::gpu_local deviceLocalKernelPrototype;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCALBASE_H_
