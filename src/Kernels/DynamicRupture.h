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

namespace seissol::kernels {

class DynamicRupture : public Kernel {
  private:
  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints m_krnlPrototype;
  kernels::Time m_timeKernel;
#ifdef ACL_DEVICE
  dynamicRupture::kernel::gpu_evaluateAndRotateQAtInterpolationPoints m_gpuKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  public:
  DynamicRupture() = default;

  void setGlobalData(const CompoundGlobalData& global) override;

  void spaceTimeInterpolation(
      const DRFaceInformation& faceInfo,
      const GlobalData* global,
      const DRGodunovData* godunovData,
      DREnergyOutput* drEnergyOutput,
      const real* timeDerivativePlus,
      const real* timeDerivativeMinus,
      real qInterpolatedPlus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
      real qInterpolatedMinus[dr::misc::TimeSteps][seissol::tensor::QInterpolated::size()],
      const real* timeDerivativePlusPrefetch,
      const real* timeDerivativeMinusPrefetch,
      const real* coeffs);

  // NOLINTNEXTLINE
  void batchedSpaceTimeInterpolation(DrConditionalPointersToRealsTable& table,
                                     const real* coeffs,
                                     seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsGodunovState(const DRFaceInformation& faceInfo,
                         std::uint64_t& nonZeroFlops,
                         std::uint64_t& hardwareFlops);
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
