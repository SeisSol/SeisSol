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
#include "Config.h"
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

class Local : public LocalKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
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
  kernel::volume<Config> volumeKernelPrototype_;
  kernel::localFlux<Config> localFluxKernelPrototype_;
  kernel::localFluxNodal<Config> nodalLfKrnlPrototype_;

  kernels::AnalyticalBoundary analyticalBoundary_;

  kernel::fsgFlux<Config> fsgFlux_;
  kernel::dirichletFlux<Config> dirichletFlux_;

#ifdef ACL_DEVICE
  kernel::gpu_volume<Config> deviceVolumeKernelPrototype_;
  kernel::gpu_localFlux<Config> deviceLocalFluxKernelPrototype_;
  kernel::gpu_localFluxAll<Config> deviceLocalFluxAllKernelPrototype_;
  kernel::gpu_localFluxNodal<Config> deviceNodalLfKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::instance();

  kernel::gpu_fsgFlux<Config> deviceFsgFlux_;
  kernel::gpu_dirichletFlux<Config> deviceDirichletFlux_;
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_
