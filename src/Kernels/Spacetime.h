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

  /**
   * @brief computeAder for the batch the recorder noted (see recording::LocalIntegrationRecorder),
   * which also applies the volume term of the corrector once it has formed the time integral: to
   * the field, or to the extended field where a solver keeps one. The device has no volume kernel,
   * so LocalKernel::computeBatchedIntegral starts from there.
   */
  virtual void computeBatchedAder(const real* coeffs,
                                  double timeStepWidth,
                                  LTS::Layer& layer,
                                  LocalTmp& tmp,
                                  recording::ConditionalPointersToRealsTable& dataTable,
                                  bool updateDisplacement,
                                  seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * @brief Write the traces of the space-time evolution on one face: what a neighbor that evaluates
   * it over intervals of its own reads of the cell there, or the fault.
   *
   * @param timeDerivativesOrSTP The space-time evolution computeAder wrote.
   * @param face The face.
   * @param side In which parametrisation of the face the traces are written
   * (CellLocalInformation::traceSides).
   * @param faceDerivatives Output: the traces, as the TimeKernel evaluates them (evaluateTrace).
   */
  virtual void computeFaceDerivatives(const real* timeDerivativesOrSTP,
                                      std::size_t face,
                                      std::size_t side,
                                      real* faceDerivatives) = 0;

  /**
   * @brief computeFaceDerivatives for the batches the recorder noted by face and by side, from the
   * space-time evolution computeBatchedAder wrote (see recording::LocalIntegrationRecorder).
   */
  virtual void
      computeBatchedFaceDerivatives(recording::ConditionalPointersToRealsTable& dataTable,
                                    seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  [[nodiscard]] virtual PerformanceEstimate metrics() const = 0;

  /// The cost of computeFaceDerivatives on one face.
  [[nodiscard]] virtual PerformanceEstimate faceDerivativesMetrics() const = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SPACETIME_H_
