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
#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"

#include <iterator>

namespace seissol::kernels::solver::linearckanelastic {
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

  /// Whether the device kernel keeps the derivatives in the kernel, for the cells that need none
  /// (cf. needsDerivatives), with the second member of its family; otherwise, all cells write
  /// their derivatives out. The code generator gives the kernel one only where it keeps the
  /// temporaries of a kernel on the chip, i.e. not with the old GPU interface. So far, it is not
  /// used (PreferKeeping): with TensorForge, the second member spills more registers than writing
  /// the derivatives costs (viscoelastic, 3 mechanisms, order 6, single precision, sm_120: 11.7
  /// instead of 8.1 KB of spill stores per thread, and 93 instead of 56 ms per step of 1e5 cells).
  static constexpr bool PreferKeeping = false;
#ifdef ACL_DEVICE
  static constexpr bool KeepsDerivativesOnDevice =
      PreferKeeping && std::size(kernel::gpu_derivative<Cfg>::ExecutePtrs) > 1;
#else
  static constexpr bool KeepsDerivativesOnDevice = false;
#endif

  [[nodiscard]] PerformanceEstimate metrics() const override;

  protected:
  kernel::derivative<Cfg> krnlPrototype_;
  kernel::fsgKernel<Cfg> fsgKernelPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_derivative<Cfg> deviceKrnlPrototype_;
#endif
};

template <typename Cfg>
class Time : public TimeKernel<Cfg> {
  public:
  using real = Real<Cfg>;

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
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIME_H_
