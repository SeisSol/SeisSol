// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_LOCAL_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_LOCAL_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Local.h"
#include "Monitoring/Metric.h"

#include <array>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::nonlinearck {

/// The cell's own contribution: the volume term and the faces it owns.
///
/// The volume term is a constant map on the transported tensor, because the
/// flux is linear in the velocity and the stress and both are transported --
/// no nodal detour, nothing to evaluate. The source integral of the internal
/// variables is added in the same kernel, being the other half of the same
/// update.
///
/// The boundary conditions are missing: they want a ghost rule for the stress
/// columns, which is where the transport layout has to reach the nodal
/// projection first.
///
/// Both differ from the linear solver in the same way. The volume term lifts a
/// flux that was evaluated at the nodes rather than applying a constant
/// operator to the state, and the faces carry a Rusanov flux rather than a
/// Riemann solution. What the two have in common is that the flux is already
/// there by the time this kernel runs: the predictor evaluated and integrated
/// it, and both the state and the stress arrive as integrals.
class Local : public LocalKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeIntegral(real* timeIntegratedDoFs,
                       LTS::Ref& data,
                       LocalTmp& tmp,
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

  [[nodiscard]] PerformanceEstimate
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const override;

  protected:
  kernel::damageCellIntegral cellIntegral_;
  kernel::projectToFace projectToFace_;
  kernel::damageRusanov rusanov_;
  kernel::faceIntegral faceIntegral_;

#ifdef ACL_DEVICE
  kernel::gpu_damageCellIntegral deviceCellIntegral_;
  kernel::gpu_projectToFace deviceProjectToFace_;
  kernel::gpu_damageRusanov deviceRusanov_;
  kernel::gpu_faceIntegral deviceFaceIntegral_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::nonlinearck

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_LOCAL_H_
