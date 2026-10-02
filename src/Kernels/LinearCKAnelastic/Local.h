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

#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Interface.h"
#include "Kernels/Local.h"
#include "Physics/InitialField.h"

#include <memory>

namespace seissol::kernels::solver::linearckanelastic {
class Local : public LocalKernel {
  public:
  void setGlobalData(const CompoundGlobalData<Config>& global) override;

  void computeIntegral(real* timeIntegratedDoFs,
                       LTS::Ref<Config>& data,
                       LocalTmp& tmp,
                       double time,
                       double timeStepWidth) override;

  void computeBatchedIntegral(recording::ConditionalPointersToRealsTable& dataTable,
                              recording::ConditionalIndicesTable& indicesTable,
                              double timeStepWidth,
                              seissol::parallel::runtime::StreamRuntime& runtime) override;

  void evaluateBatchedTimeDependentBc(recording::ConditionalPointersToRealsTable& dataTable,
                                      recording::ConditionalIndicesTable& indicesTable,
                                      LTS::Layer& layer,
                                      double time,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] PerformanceEstimate
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const override;

  protected:
  kernel::volumeExt<Config> volumeKernelPrototype_;
  kernel::localFluxExt<Config> localFluxKernelPrototype_;
  kernel::local<Config> localKernelPrototype_;

  kernel::fsgFlux<Config> fsgFlux_;
  kernel::dirichletFlux<Config> dirichletFlux_;
  kernel::localFluxNodal<Config> nodalLfKrnlPrototype_;

  kernels::AnalyticalBoundary analyticalBoundary_;

#ifdef ACL_DEVICE
  kernel::gpu_volumeExt<Config> deviceVolumeKernelPrototype_;
  kernel::gpu_localFluxExt<Config> deviceLocalFluxKernelPrototype_;
  kernel::gpu_local<Config> deviceLocalKernelPrototype_;
  kernel::gpu_fluxLocalAll<Config> deviceFluxLocalAllKernelPrototype_;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCAL_H_
