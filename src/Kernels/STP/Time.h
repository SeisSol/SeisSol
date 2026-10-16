// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_STP_TIME_H_
#define SEISSOL_SRC_KERNELS_STP_TIME_H_

#include "Common/Constants.h"
#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"
#include "Monitoring/Metric.h"

#include <iterator>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol::kernels::solver::stp {

template <typename Cfg>
class Spacetime : public SpacetimeKernel<Cfg> {
  public:
  using real = Real<Cfg>;

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

  /// Whether the device kernel keeps the predictor in the kernel, for the cells that need none
  /// (cf. needsDerivatives), with the second member of its family; otherwise, all cells write
  /// their predictor out. The code generator gives the kernel one only where it keeps the
  /// temporaries of a kernel on the chip, i.e. not with the old GPU interface. So far, it is not
  /// used (PreferKeeping): with TensorForge, keeping the predictor on the chip costs more than
  /// writing it out (poroelastic, order 6, single precision, sm_120: 155 instead of 111 ms per step
  /// of 1e5 cells).
  static constexpr bool PreferKeeping = false;
#ifdef ACL_DEVICE
  static constexpr bool KeepsDerivativesOnDevice =
      PreferKeeping && std::size(kernel::gpu_spaceTimePredictor<Cfg>::ExecutePtrs) > 1;
#else
  static constexpr bool KeepsDerivativesOnDevice = false;
#endif

  [[nodiscard]] PerformanceEstimate metrics() const override;

  private:
  void executeSTP(double timeStepWidth, LTS::Ref<Cfg>& data, real* timeIntegrated, real* stp);

  kernel::spaceTimePredictor<Cfg> krnlPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_spaceTimePredictor<Cfg> deviceKrnlPrototype_;
#endif
};

template <typename Cfg>
class Time : public TimeKernel<Cfg> {
  public:
  using real = Real<Cfg>;

  void setGlobalData(const CompoundGlobalData<Cfg>& global) override;
  void evaluate(const real* coeffs, const real* timeDerivatives, real* timeEvaluated) override;
  void evaluateBatched(const real* coeffs,
                       const real** timeDerivatives,
                       real** timeIntegratedDofs,
                       std::size_t numElements,
                       seissol::parallel::runtime::StreamRuntime& runtime) override;
  [[nodiscard]] PerformanceEstimate metrics() const override;
};

} // namespace seissol::kernels::solver::stp

#endif // SEISSOL_SRC_KERNELS_STP_TIME_H_
