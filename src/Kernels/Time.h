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
#include <cstddef>
#include <memory>

namespace seissol::kernels {

class TimeKernel : public Kernel {
  public:
  ~TimeKernel() override = default;

  /**
    @brief Evaluates a given space-time representation in time.

    Can be used for both point evaluation, derivative computation, and integration.

    Is best used in conjunction with the results from a TimeBasis class function.

    @param coeffs The time basis coefficients.
    @param timeDerivatives A pointer to the input space-time data.
    @param timeEvaluated A pointer to the returned time-evaluated data.
  */
  virtual void evaluate(const real* coeffs,
                        const real* timeDerivatives,
                        real timeEvaluated[tensor::I::size()]) = 0;

  /**
    @brief Evaluates a given space-time representation in time.

    Cf. the comments to `evaluate`.

    @param coeffs The time basis coefficients.
    @param timeDerivatives A batch pointer to the input space-time data.
    @param timeEvaluated A batch pointer to the returned time-evaluated data.
    @param numElements The number of elements to process (indicates the number of elements in
    timeDerivatives and timeEvaluated)
    @param runtime The stream to place the operations on.
  */
  virtual void evaluateBatched(const real* coeffs,
                               const real** timeDerivatives,
                               real** timeIntegratedDofs,
                               std::size_t numElements,
                               seissol::parallel::runtime::StreamRuntime& runtime) = 0;
  virtual void flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIME_H_
