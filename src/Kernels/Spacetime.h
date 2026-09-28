// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_SPACETIME_H_
#define SEISSOL_SRC_KERNELS_SPACETIME_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/Kernel.h"
#include "Monitoring/Metric.h"
#include "Numerical/BasisFunction.h"
#include "Parallel/Runtime/Stream.h"

#include <cassert>
#include <limits>
#include <memory>

namespace seissol::kernels {

class SpacetimeKernel : public Kernel {
  public:
  ~SpacetimeKernel() override = default;

  /**
   * @brief Compute the space-time evolution (a.k.a. "derivatives") from the given DOFs.
   *
   * Includes an inline time integration; as that's needed for most subsequent operations.
   *
   * This step in essence equals the ADER-DG "predictor" step. Note that the exact method
   * for computing the space-time evolution may vary (e.g. CK vs STP); see the respective subfolder
   * names for more infos.
   *
   * @param coeffs The time basis coefficients used for the integrated integration.
   * @param timeStepWidth The size of the current timestep
   * @param data Cell data reference object (contains references to all stored data arrays for that
   * cell)
   * @param tmp Local integration temporary data object (contains thread-local data for a cell; to
   * pass to the integration computation)
   * @param timeIntegrated Output: time integration data.
   * @param timeDerivativesOrSTP Output: space-time evoluion.
   * @param updateDisplacement Update the face displacement (needed for elastic-acoustic)
   */
  virtual void computeAder(const real* coeffs,
                           double timeStepWidth,
                           LTS::Ref& data,
                           LocalTmp& tmp,
                           real* timeIntegrated,
                           real* timeDerivativesOrSTP = nullptr,
                           bool updateDisplacement = false) = 0;

  virtual void computeBatchedAder(const real* coeffs,
                                  double timeStepWidth,
                                  LTS::Layer& layer,
                                  LocalTmp& tmp,
                                  recording::ConditionalPointersToRealsTable& dataTable,
                                  recording::ConditionalMaterialTable& materialTable,
                                  bool updateDisplacement,
                                  seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * @brief Enqueue the host-side preparations of computeBatchedAder (if any).
   *
   * Needs to be called before each computeBatchedAder call with the same time step width, and
   * outside of any device graph capture (since host work cannot be captured on all backends;
   * e.g. on HIP, enqueueHost synchronizes the stream).
   *
   * @param timeStepWidth The size of the current timestep
   * @param layer The layer which computeBatchedAder is called on
   * @param runtime The stream runtime to enqueue the preparations on
   */
  virtual void prepareBatchedAder(double timeStepWidth,
                                  LTS::Layer& layer,
                                  seissol::parallel::runtime::StreamRuntime& runtime) {}

  [[nodiscard]] virtual PerformanceEstimate metrics() const = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SPACETIME_H_
