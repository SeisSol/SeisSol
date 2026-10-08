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

#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/AnalyticalBoundary.h"
#include "Kernels/Interface.h"
#include "Kernels/Local.h"
#include "Physics/InitialField.h"

#include <memory>

namespace seissol::kernels::solver::linearckanelastic {
template <typename Cfg>
class Local : public LocalKernel<Cfg> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  void setGlobalData(const CompoundGlobalData<Cfg>& global) override;

  void computeIntegral(real* timeIntegratedDoFs,
                       LTS::Ref<Cfg>& data,
                       LocalTmp<Cfg>& tmp,
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
  kernel::volumeExt<Cfg> volumeKernelPrototype_;
  kernel::localFluxExt<Cfg> localFluxKernelPrototype_;
  kernel::local<Cfg> localKernelPrototype_;

  kernel::fsgFlux<Cfg> fsgFlux_;
  kernel::dirichletFlux<Cfg> dirichletFlux_;
  kernel::localFluxNodal<Cfg> nodalLfKrnlPrototype_;

  kernels::AnalyticalBoundary<Cfg> analyticalBoundary_;

#ifdef ACL_DEVICE
  kernel::gpu_volumeExt<Cfg> deviceVolumeKernelPrototype_;
  kernel::gpu_localFluxExt<Cfg> deviceLocalFluxKernelPrototype_;
  kernel::gpu_local<Cfg> deviceLocalKernelPrototype_;
  kernel::gpu_fluxLocalAll<Cfg> deviceFluxLocalAllKernelPrototype_;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_LOCAL_H_
