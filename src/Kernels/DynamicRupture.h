// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_

#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Kernel.h"
#include "Kernels/Solver.h"
#include "Monitoring/Metric.h"

namespace seissol::kernels {

/// The state of both sides of a fault face at its points and the points in time the friction law
/// is evaluated at, from the traces of their time derivatives on the face (LTS::FaceBuffers). Both
/// sides write those in the parametrisation of the plus side, so one kernel reads either side and
/// any face of a cell.
class DynamicRupture : public Kernel {
  private:
  dynamicRupture::kernel::evaluateAndRotateTraceAtInterpolationPoints krnlPrototype_;
  kernels::Time timeKernel_;
#ifdef ACL_DEVICE
  dynamicRupture::kernel::gpu_projectTracesToDR gpuCombinedKrnlPrototype_;
#endif

  public:
  DynamicRupture() = default;

  void setGlobalData(const CompoundGlobalData& global) override;

  /// @param timeDerivativePlus the traces of the time derivatives of the plus side on the face
  /// @param timeDerivativeMinus the same for the minus side
  void spaceTimeInterpolation(
      const DRGodunovData* godunovData,
      const real* timeDerivativePlus,
      const real* timeDerivativeMinus,
      real qInterpolatedPlus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
      real qInterpolatedMinus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
      const real* coeffs);

  // NOLINTNEXTLINE
  void batchedSpaceTimeInterpolation(recording::DrConditionalPointersToRealsTable& table,
                                     const real* coeffs,
                                     seissol::parallel::runtime::StreamRuntime& runtime);

  /// what one fault face takes
  [[nodiscard]] PerformanceEstimate metrics() const;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
