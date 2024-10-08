/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Dynamic Rupture kernel of SeisSol.
 **/

#ifndef KERNELS_DYNAMICRUPTURE_H_
#define KERNELS_DYNAMICRUPTURE_H_

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

#endif
