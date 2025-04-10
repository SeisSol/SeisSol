// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_TIME_H_
#define SEISSOL_SRC_KERNELS_TIME_H_

#include "Initializer/Typedefs.h"
#include "generated_code/tensor.h"
#include <Kernels/Kernel.h>
#include <Numerical/BasisFunction.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <memory>

namespace seissol::kernels {

class TimeKernel : public Kernel {
  public:
  ~TimeKernel() override = default;

  virtual void evaluateAtTime(
      std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions,
      const real* timeDerivatives,
      real timeEvaluated[tensor::Q::size()]) = 0;
  virtual void flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops) = 0;

  virtual void computeIntegral(double expansionPoint,
                               double integrationStart,
                               double integrationEnd,
                               const real* timeDerivatives,
                               real timeIntegrated[tensor::I::size()]) = 0;

  virtual void computeBatchedIntegral(double expansionPoint,
                                      double integrationStart,
                                      double integrationEnd,
                                      const real** timeDerivatives,
                                      real** timeIntegratedDofs,
                                      unsigned numElements,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void computeTaylorExpansion(real time,
                                      real expansionPoint,
                                      const real* timeDerivatives,
                                      real timeEvaluated[tensor::Q::size()]) = 0;

  virtual void
      computeBatchedTaylorExpansion(real time,
                                    real expansionPoint,
                                    real** timeDerivatives,
                                    real** timeEvaluated,
                                    size_t numElements,
                                    seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIME_H_
