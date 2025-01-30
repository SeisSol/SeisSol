// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Time.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#ifdef USE_STP
#include "Numerical/BasisFunction.h"
#include <array>
#include <memory>
#endif

namespace seissol::kernels {

constexpr std::size_t NumSpaceQuadraturePoints = (ConvergenceOrder + 1) * (ConvergenceOrder + 1);

class DynamicRupture {
  private:
  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints m_krnlPrototype;
  kernels::Time m_timeKernel;
#ifdef ACL_DEVICE
  dynamicRupture::kernel::gpu_evaluateAndRotateQAtInterpolationPoints m_gpuKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  public:
  double timePoints[ConvergenceOrder]{};
  double timeWeights[ConvergenceOrder]{};
  real spaceWeights[NumSpaceQuadraturePoints]{};
#ifdef USE_STP
  std::array<std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>>, ConvergenceOrder>
      timeBasisFunctions;
#endif

  DynamicRupture() = default;

  static void checkGlobalData(const GlobalData* global, size_t alignment);
  void setHostGlobalData(const GlobalData* global);
  void setGlobalData(const CompoundGlobalData& global);

  void setTimeStepWidth(double timestep);

  void spaceTimeInterpolation(
      const DRFaceInformation& faceInfo,
      const GlobalData* global,
      const DRGodunovData* godunovData,
      DREnergyOutput* drEnergyOutput,
      const real* timeDerivativePlus,
      const real* timeDerivativeMinus,
      real qInterpolatedPlus[ConvergenceOrder][seissol::tensor::QInterpolated::size()],
      real qInterpolatedMinus[ConvergenceOrder][seissol::tensor::QInterpolated::size()],
      const real* timeDerivativePlusPrefetch,
      const real* timeDerivativeMinusPrefetch);

  // NOLINTNEXTLINE
  void batchedSpaceTimeInterpolation(DrConditionalPointersToRealsTable& table,
                                     seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsGodunovState(const DRFaceInformation& faceInfo,
                         long long& nonZeroFlops,
                         long long& hardwareFlops);
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_DYNAMICRUPTURE_H_
