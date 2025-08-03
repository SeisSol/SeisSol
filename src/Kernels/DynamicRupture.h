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
#include <Kernels/Kernel.h>
#include <Kernels/Solver.h>

namespace seissol::kernels {

template<typename Cfg>
class DynamicRupture : public Kernel {
  private:
  using real = Real<Cfg>;

  dynamicRupture::kernel::evaluateAndRotateQAtInterpolationPoints<Cfg> m_krnlPrototype;
  kernels::Time<Cfg> m_timeKernel;
#ifdef ACL_DEVICE
  dynamicRupture::kernel::gpu_evaluateAndRotateQAtInterpolationPoints<Cfg> m_gpuKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  public:
  DynamicRupture() = default;

  void setGlobalData(const CompoundGlobalData& global) override;

  void spaceTimeInterpolation(
      const DRFaceInformation& faceInfo,
      const GlobalData* global,
      const DRGodunovData<Cfg>* godunovData,
      DREnergyOutput<Cfg>* drEnergyOutput,
      const real* timeDerivativePlus,
      const real* timeDerivativeMinus,
      real qInterpolatedPlus[Cfg::ConvergenceOrder][seissol::tensor::QInterpolated<Cfg>::size()],
      real qInterpolatedMinus[Cfg::ConvergenceOrder][seissol::tensor::QInterpolated<Cfg>::size()],
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
