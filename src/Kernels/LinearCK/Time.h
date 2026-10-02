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
#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"
#include "Monitoring/Metric.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
class Spacetime : public SpacetimeKernel<Cfg> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  void setGlobalData(const CompoundGlobalData<Cfg>& global) override;
  void computeAder(const real* coeffs,
                   double timeStepWidth,
                   LTS::Ref<Cfg>& data,
                   LocalTmp<Cfg>& tmp,
                   real* timeIntegrated,
                   real* timeDerivativesOrSTP = nullptr,
                   bool updateDisplacement = false) override;
  void computeBatchedAder(const real* coeffs,
                          double timeStepWidth,
                          LTS::Layer& layer,
                          LocalTmp<Cfg>& tmp,
                          recording::ConditionalPointersToRealsTable& dataTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] PerformanceEstimate metrics() const override;

  protected:
  kernel::derivative<Cfg> krnlPrototype_;

  kernel::fsgKernel<Cfg> fsgKernelPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_derivative<Cfg> deviceKrnlPrototype_;
  kernel::gpu_fsgKernel<Cfg> deviceFsgKernelPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::instance();
#endif
};

template <typename Cfg>
class Time : public TimeKernel<Cfg> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  void setGlobalData(const CompoundGlobalData<Cfg>& global) override;
  void evaluate(const real* coeffs,
                const real* timeDerivatives,
                real timeEvaluated[tensor::I<Cfg>::size()]) override;
  void evaluateBatched(const real* coeffs,
                       const real** timeDerivatives,
                       real** timeIntegratedDofs,
                       std::size_t numElements,
                       seissol::parallel::runtime::StreamRuntime& runtime) override;
  [[nodiscard]] PerformanceEstimate metrics() const override;
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_
