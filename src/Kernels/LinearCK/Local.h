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
#include "Kernels/CurvedSplit.h"
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
                       LTS::Ref& data,
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

  [[nodiscard]] PerformanceEstimate
      curvedMetrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const override {
    return decltype(curvedSplit_)::metrics(faceTypes);
  }

  protected:
  kernel::volume volumeKernelPrototype_;
  kernel::faceTrace faceTraceKernelPrototype_;
  kernel::localFluxNodal nodalLfKrnlPrototype_;
  // the volume term of a curved cell and what its momentum rows take back on its faces
  kernels::CurvedSplitT curvedSplit_;

  kernels::AnalyticalBoundary analyticalBoundary_;

  kernel::fsgFlux fsgFlux_;
  kernel::dirichletFlux dirichletFlux_;
  // what a boundary folds into the local flux solver, where the corrector does not apply it
  // (appliesFoldApart)
  kernel::foldedBoundaryFlux foldedBoundaryFlux_;

#ifdef ACL_DEVICE
  // the traces of the four faces of a cell at once, with the side of each face as data
  kernel::gpu_faceTraces deviceFaceTracesKernelPrototype_;
  kernel::gpu_localFluxNodal deviceNodalLfKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();

  kernel::gpu_fsgFlux deviceFsgFlux_;
  kernel::gpu_dirichletFlux deviceDirichletFlux_;
  kernel::gpu_foldedBoundaryFlux deviceFoldedBoundaryFlux_;
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_LOCAL_H_
