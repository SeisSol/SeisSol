// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_TIME_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_TIME_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"
#include "Monitoring/Metric.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::nonlinearck {

/// The predictor for a material whose flux is nonlinear in the state.
///
/// The Cauchy-Kovalevskaya recursion runs unchanged, on star matrices that
/// linearise the material about the cell mean. What it produces is a time
/// expansion of the state, not of the flux -- for a nonlinear flux the two are
/// not the same thing. The flux is therefore evaluated pointwise at the time
/// quadrature nodes and integrated there, which is also where the internal
/// variables pick up their source terms.
///
/// Two integrals leave this kernel: the state and the stress. The stress is
/// carried because the integral of a nonlinear flux is not the flux of the
/// integrated state, and because carrying it spares every neighbour the
/// material of the cell it reads from.
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
                          recording::ConditionalMaterialTable& materialTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] PerformanceEstimate metrics() const override;

  protected:
  kernel::derivative derivative_;
  kernel::convertToNodal convertToNodal_;
  kernel::damageInvariants invariants_;
  kernel::damageStress stress_;
  kernel::damageFlux flux_;
  kernel::damageCellState cellState_;
  kernel::damageSource source_;
  kernel::projectDerivativeToNodalBoundaryRotated projectDerivativeToNodalBoundaryRotated_;

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceDerivative_;
  kernel::gpu_convertToNodal deviceConvertToNodal_;
  kernel::gpu_damageInvariants deviceInvariants_;
  kernel::gpu_damageStress deviceStress_;
  kernel::gpu_damageFlux deviceFlux_;
  kernel::gpu_damageCellState deviceCellState_;
  kernel::gpu_damageSource deviceSource_;
  kernel::gpu_projectDerivativeToNodalBoundaryRotated deviceDerivativeToNodalBoundaryRotated_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif
};

/// Evaluation of the stored expansion at a point in time. The expansion is of
/// the state, so this is the same Taylor sum the linear solver uses.
class Time : public TimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void evaluate(const real* coeffs,
                const real* timeDerivatives,
                real timeEvaluated[tensor::I::size()]) override;
  void evaluateBatched(const real* coeffs,
                       const real** timeDerivatives,
                       real** timeIntegratedDofs,
                       std::size_t numElements,
                       seissol::parallel::runtime::StreamRuntime& runtime) override;
  [[nodiscard]] PerformanceEstimate metrics() const override;
};

} // namespace seissol::kernels::solver::nonlinearck

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_TIME_H_
