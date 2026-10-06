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
#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Local.h"
#include "Monitoring/Metric.h"

#include <memory>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "Kernels/AnalyticalBoundary.h"
#pragma GCC diagnostic pop
#include "Physics/InitialField.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
class Local : public LocalKernel<Cfg> {
  public:
  using real = Real<Cfg>;

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
  kernel::volume<Cfg> volumeKernelPrototype_;
  kernel::localFlux<Cfg> localFluxKernelPrototype_;
  kernel::localFluxNodal<Cfg> nodalLfKrnlPrototype_;

  kernels::AnalyticalBoundary<Cfg> analyticalBoundary_;

  kernel::fsgFlux<Cfg> fsgFlux_;
  kernel::dirichletFlux<Cfg> dirichletFlux_;

#ifdef ACL_DEVICE
  kernel::gpu_volume<Cfg> deviceVolumeKernelPrototype_;
  kernel::gpu_localFlux<Cfg> deviceLocalFluxKernelPrototype_;
  kernel::gpu_localFluxAll<Cfg> deviceLocalFluxAllKernelPrototype_;
  kernel::gpu_localFluxNodal<Cfg> deviceNodalLfKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::instance();

  kernel::gpu_fsgFlux<Cfg> deviceFsgFlux_;
  kernel::gpu_dirichletFlux<Cfg> deviceDirichletFlux_;
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_
