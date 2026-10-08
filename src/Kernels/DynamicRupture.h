// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_

#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Kernel.h"
#include "Kernels/Solver.h"
#include "Monitoring/Metric.h"

namespace seissol::kernels {

/// The interpolation of the space-time predictor onto the faults, for the configuration `Cfg`.
template <typename Cfg>
class DynamicRupture : public Kernel<Cfg> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  private:
  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg> krnlPrototype_;
  kernels::Time<Cfg> timeKernel_;
#ifdef ACL_DEVICE
  dynamicRupture::kernel::gpu_evaluateAndRotateQAtInterpolationPoints<Cfg> gpuKrnlPrototype_;
  dynamicRupture::kernel::gpu_projectToDR<Cfg> gpuCombinedKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::instance();
#endif

  public:
  DynamicRupture() = default;

  void setGlobalData(const CompoundGlobalData<Cfg>& global) override;

  void spaceTimeInterpolation(
      const DRFaceInformation& faceInfo,
      const DRGodunovData<Cfg>* godunovData,
      const real* timeDerivativePlus,
      const real* timeDerivativeMinus,
      real qInterpolatedPlus[dr::misc::TimeSteps<Cfg>][seissol::tensor::QInterpolated<Cfg>::size()],
      real qInterpolatedMinus[dr::misc::TimeSteps<Cfg>]
                             [seissol::tensor::QInterpolated<Cfg>::size()],
      const real* timeDerivativePlusPrefetch,
      const real* timeDerivativeMinusPrefetch,
      const real* coeffs);

  // NOLINTNEXTLINE
  void batchedSpaceTimeInterpolation(recording::DrConditionalPointersToRealsTable& table,
                                     const real* coeffs,
                                     seissol::parallel::runtime::StreamRuntime& runtime);

  [[nodiscard]] PerformanceEstimate metrics(const DRFaceInformation& faceInfo) const;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
