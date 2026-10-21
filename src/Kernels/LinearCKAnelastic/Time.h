// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIME_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIME_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"

namespace seissol::kernels::solver::linearckanelastic {
class Spacetime : public SpacetimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeAder(const real* coeffs,
                   double timeStepWidth,
                   LTS::Ref& data,
                   LocalTmp& tmp,
                   real* timeIntegrated,
                   real* timeDerivativesOrSTP = nullptr,
                   bool updateDisplacement = false) override;
  void computeBatchedAder(const real* coeffs,
                          double timeStepWidth,
                          LTS::Layer& layer,
                          LocalTmp& tmp,
                          recording::ConditionalPointersToRealsTable& dataTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime) override;

  void computeFaceDerivatives(const real* timeDerivativesOrSTP,
                              std::size_t face,
                              std::size_t side,
                              real* faceDerivatives) override;
  void computeBatchedFaceDerivatives(recording::ConditionalPointersToRealsTable& dataTable,
                                     seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] PerformanceEstimate metrics() const override;
  [[nodiscard]] PerformanceEstimate faceDerivativesMetrics() const override;

  protected:
  kernel::derivative krnlPrototype_;
  kernel::faceDerivatives faceDerivativesKrnlPrototype_;
  kernel::fsgKernel fsgKernelPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceKrnlPrototype_;
  kernel::gpu_faceDerivatives deviceFaceDerivativesKrnlPrototype_;
#endif
};

class Time : public TimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void evaluate(const real* coeffs,
                const real* timeDerivatives,
                real timeEvaluated[tensor::I::size()]) override;
  void evaluateTrace(const real* coeffs, const real* faceDerivatives, real* faceTrace) override;
  void evaluateTraceBatched(const real* coeffs,
                            const real** faceDerivatives,
                            real** faceTraces,
                            std::size_t numElements,
                            seissol::parallel::runtime::StreamRuntime& runtime) override;
  [[nodiscard]] PerformanceEstimate metrics() const override;
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIME_H_
