// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Time.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/NonLinearCK/Solver.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Numerical/TimeBasis.h"
#include "Parallel/Runtime/Stream.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utils/logger.h>
#include <yateto.h>
#include <yateto/InitTools.h>

namespace seissol::kernels::solver::nonlinearck {

namespace {
/// Below this, the second strain invariant is taken for zero. It guards a
/// square root and a division, so it is tied to the working precision and not
/// to the material: the square of the machine epsilon is far below any strain
/// a simulation resolves and far above where a reciprocal stops being finite.
constexpr real invariantFloor() {
  return std::numeric_limits<real>::epsilon() * std::numeric_limits<real>::epsilon();
}
} // namespace

void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  derivative_.bindGlobals(*global.onHost);
  step_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceDerivative_.bindGlobals(*global.onDevice);
  deviceStep_.bindGlobals(*global.onDevice);
#endif
}

void Spacetime::computeAder(const real* coeffs,
                            double timeStepWidth,
                            LTS::Ref& data,
                            LocalTmp& tmp,
                            real* timeIntegrated,
                            real* timeDerivatives,
                            bool /*updateDisplacement*/) {
  assert(reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>()) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(timeIntegrated) % Alignment == 0);
  assert(timeDerivatives == nullptr ||
         reinterpret_cast<uintptr_t>(timeDerivatives) % Alignment == 0);

  alignas(PagesizeStack) real temporaryBuffer[Solver::DerivativesSize];
  auto* derivativesBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;

  const auto& local = data.get<LTS::LocalIntegration>().specific;

  // The expansion of the state, and with it the columns of the transported
  // tensor that the state couples through.
  kernel::derivative derivative = derivative_;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    derivative.star(i) = data.get<LTS::LocalIntegration>().starMatrices[i];
  }
  derivative.dQ(0) = const_cast<real*>(data.get<LTS::Dofs>());
  for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    derivative.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  derivative.I = timeIntegrated;
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    derivative.power(der) = coeffs[der];
  }
  derivative.execute();

  // Everything nonlinear, in one kernel. The rule it samples with is ours to
  // choose: the nodes and weights of the quadrature, the coefficients that
  // evaluate the expansion there, and how far the internal variables march
  // from one node to the next.
  const Solver::TimeBasis<real> basis(ConvergenceOrder);
  const auto [nodes, weights] = basis.quadratureWithEndpoints(timeStepWidth);

  kernel::damageStep step = step_;
  step.dQ(0) = data.get<LTS::Dofs>();
  for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    step.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  for (std::size_t q = 0; q < nodes.size(); ++q) {
    const auto evaluation = basis.point(nodes[q], timeStepWidth);
    for (std::size_t i = 0; i < ConvergenceOrder; ++i) {
      step.evaluate(q, i) = evaluation[i];
    }
    step.weight(q) = weights[q];
    // The internal variables are carried across the nodes explicitly. The
    // first node is the start of the step and the last one its end, so the
    // marches tile the step without a gap; what is left is the order of the
    // march itself, and that is a property of these numbers alone.
    step.march(q) = (q + 1 < nodes.size() ? nodes[q + 1] : nodes[q]) - nodes[q];
  }

  step.I = timeIntegrated;
  step.sourceI = tmp.sourceIntegral;
  step.epsInit = local.epsInit;

  step.materialParameters = local.parameters;
  step.invariantFloor = invariantFloor();
  step.execute();

  // What scales the dissipation is in the transported tensor as well, and
  // the step wrote it there itself: the square of the fastest wave integrated
  // over the step, and the length of the step beside it. Neither is a speed,
  // which is what lets a coarser cluster accumulate them like every other
  // column and a face divide one by the other.
}

void Spacetime::computeBatchedAder(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM LocalTmp& tmp,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalMaterialTable& materialTable,
    SEISSOL_GPU_PARAM bool updateDisplacement,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
  // Everything this needs is recorded: the expansion, the transported tensor,
  // the source integral, and the initial strain through an offset into the
  // cell's local integration data. What is not is the material: the step
  // kernel takes it as scalars, and a scalar argument is uniform over a
  // batch, so the thirteen of them have to be read per element the way the
  // face's wave speed already is.
  logError() << "The batched step kernel needs the material per element, not as scalars.";
}

PerformanceEstimate Spacetime::metrics() const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::derivative>();
  estimate += PerformanceEstimate::fromKernel<kernel::damageStep>();

  std::uint64_t reals = 0;
  // the state in, the transported tensor and the source integral out
  reals += tensor::Q::size() + tensor::I::size() + tensor::sourceI::size();
  // the star matrices and the expansion the step reads back
  reals += yateto::computeFamilySize<tensor::star>();
  reals += yateto::computeFamilySize<tensor::dQ>();

  estimate.bytes = reals * sizeof(real);
  return estimate;
}

void Time::evaluate(const real* coeffs,
                    const real* timeDerivatives,
                    real timeEvaluated[tensor::I::size()]) {
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Alignment == 0);

  // The expansion is of the state, so this is the Taylor sum the linear
  // solver evaluates -- over the columns the two tensors share. The stress
  // has no expansion stored, so a subinterval of it cannot be reconstructed
  // here; that is what SupportsLTS being false says.
  kernel::derivativeTaylorExpansion krnl;
  krnl.I = timeEvaluated;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.execute();
}

void Time::evaluateBatched(SEISSOL_GPU_PARAM const real* coeffs,
                           SEISSOL_GPU_PARAM const real** timeDerivatives,
                           SEISSOL_GPU_PARAM real** timeIntegratedDofs,
                           SEISSOL_GPU_PARAM std::size_t numElements,
                           SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

  // The expansion is of the state, so this is the Taylor sum of the linear
  // solver, over the columns the two tensors share. What it does not write is
  // the stress, which has no expansion stored, and that is what SupportsLTS
  // being false says.
  kernel::gpu_derivativeTaylorExpansion krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = timeDerivatives;
    krnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate Time::metrics() const {
  return PerformanceEstimate::fromKernel<kernel::derivativeTaylorExpansion>();
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::nonlinearck
