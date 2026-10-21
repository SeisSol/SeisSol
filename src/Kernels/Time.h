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
#include "Kernels/Kernel.h"
#include "Monitoring/Metric.h"
#include "Numerical/BasisFunction.h"
#include "Parallel/Runtime/Stream.h"

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
  virtual void evaluate(const real* coeffs, const real* timeDerivatives, real* timeEvaluated) = 0;

  /**
    @brief Evaluates the traces of a space-time representation on a face in time, as `evaluate`
    does in the volume, with the same coefficients.

    @param coeffs The time basis coefficients.
    @param faceDerivatives A pointer to the traces (SpacetimeKernel::computeFaceDerivatives).
    @param faceTrace A pointer to the returned time-evaluated trace.
  */
  virtual void evaluateTrace(const real* coeffs, const real* faceDerivatives, real* faceTrace) = 0;

  /**
    @brief Evaluates the traces of space-time representations on faces in time.

    Cf. the comments to `evaluateTrace`.

    @param coeffs The time basis coefficients.
    @param faceDerivatives A batch pointer to the traces.
    @param faceTraces A batch pointer to the returned time-evaluated traces.
    @param numElements The number of traces to evaluate.
    @param runtime The stream to place the operations on.
  */
  virtual void evaluateTraceBatched(const real* coeffs,
                                    const real** faceDerivatives,
                                    real** faceTraces,
                                    std::size_t numElements,
                                    seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  [[nodiscard]] virtual PerformanceEstimate metrics() const = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIME_H_
