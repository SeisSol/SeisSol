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

#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
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

  virtual void evaluate(const real* coeffs,
                        const real* timeDerivatives,
                        real timeEvaluated[tensor::I::size()]) = 0;
  virtual void evaluateBatched(const real* coeffs,
                               const real** timeDerivatives,
                               real** timeIntegratedDofs,
                               std::size_t numElements,
                               seissol::parallel::runtime::StreamRuntime& runtime) = 0;
  virtual void flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIME_H_
