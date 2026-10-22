// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"
#include "Monitoring/Metric.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::kernels::solver::linearck {

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
  kernel::gpu_fsgKernel deviceFsgKernelPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
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

  private:
  // the Taylor sum of the derivatives, with the global matrices it reads bound (the elevations
  // of the degrees where the derivatives are kept in the bases of their degrees)
  kernel::derivativeTaylorExpansion evaluateKrnlPrototype_;
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_
