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
#include "Common/Offset.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/NonLinearCK/Solver.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
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

/// How many time nodes the step kernel was generated for. Its table of
/// evaluation coefficients is indexed (node, coefficient), so the stride from
/// one coefficient to the next is the number of nodes.
constexpr std::size_t stepNodes() { return tensor::evaluate::index(0, 1); }
static_assert(stepNodes() == ConvergenceOrder + 1,
              "The step samples with a rule that has the ends of the step among its nodes, which "
              "is one node more than the order. The generated kernel was built for a different "
              "number, and writing the rule into it would write past its tables.");
} // namespace

void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  derivative_.bindGlobals(*global.onHost);
  transport_.bindGlobals(*global.onHost);
  step_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceDerivative_.bindGlobals(*global.onDevice);
  deviceTransport_.bindGlobals(*global.onDevice);
  deviceStep_.bindGlobals(*global.onDevice);
#endif
}

void Spacetime::computeAder(const TimeStepCoefficients& coeffs,
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

  // The recursion reads the state where the cell keeps it, so the slot the
  // expansion reserves for its zeroth coefficient stays empty here. Whoever
  // reads the expansion back reaches that slot like any other, and would take
  // whatever the buffer happened to hold for the state -- so the state goes
  // there, as the linear solver puts it there for the same reason.
  if (timeDerivatives != nullptr) {
    std::copy_n(data.get<LTS::Dofs>(), tensor::dQ::size(0), derivativesBuffer);
  }

  const auto& local = data.get<LTS::LocalIntegration>().specific;

  // What the recursion transports by, at the cell mean of this step. The
  // geometry it is built from is what the cell keeps.
  kernel::damageTransport assemble = transport_;
  assemble.Q = data.get<LTS::Dofs>();
  assemble.epsInit = local.epsInit;
  assemble.materialParameters = local.parameters;
  assemble.invariantFloor = invariantFloor();
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    assemble.star(i) = data.get<LTS::LocalIntegration>().starMatrices[i];
    assemble.transport(i) = tmp.transport[i];
  }
  assemble.execute();

  // The expansion of the state, and with it the columns of the transported
  // tensor that the state couples through.
  kernel::derivative derivative = derivative_;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::transport>(); ++i) {
    derivative.transport(i) = tmp.transport[i];
  }
  derivative.dQ(0) = const_cast<real*>(data.get<LTS::Dofs>());
  for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    derivative.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  derivative.I = timeIntegrated;
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    derivative.power(der) = coeffs.integral.state[der];
  }
  derivative.execute();

  // Everything nonlinear, in one kernel. Where it samples is a question about
  // the step rather than about this cell, so the rule arrives with the step:
  // the nodes, their weights, and the coefficients that evaluate the
  // expansions at each of them.
  const auto& quadrature = coeffs.quadrature;

  kernel::damageStep step = step_;
  step.dQ(0) = data.get<LTS::Dofs>();
  for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    step.dQ(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  assert(quadrature.nodes.size() == stepNodes());
  for (std::size_t q = 0; q < quadrature.nodes.size(); ++q) {
    for (std::size_t i = 0; i < ConvergenceOrder; ++i) {
      step.evaluate(q, i) = quadrature.coefficients[q].state[i];
    }
    step.weight(q) = quadrature.weights[q];
  }
  // How far the internal variables have traveled by each node follows from
  // where the nodes are, which the kernel knows; how wide the step those
  // nodes span is, it does not.
  step.stepWidth = timeStepWidth;

  step.I = timeIntegrated;
  step.sourceI = tmp.sourceIntegral;
  step.epsInit = local.epsInit;

  // The expansion of what the cell carries beyond its state sits behind the
  // expansion of the state, in the same buffer: both are what a neighbour on
  // a coarser cluster reads, and both are written here.
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::transportDer>(); ++i) {
    step.transportDer(i) = derivativesBuffer + yateto::computeFamilySize<tensor::dQ>() +
                           yateto::computeFamilySize<tensor::transportDer>(1, i);
  }

  step.materialParameters = local.parameters;
  step.invariantFloor = invariantFloor();
  step.execute();

  // What scales the dissipation is in the transported tensor as well, and the
  // step wrote it there itself: the largest square of each wave family's speed
  // over the step. The squares, because that is what the moduli are affine in,
  // and the largest rather than a mean because Rusanov wants a bound -- which
  // is also what lets a coarser cluster fold one step into another by taking
  // the larger of the two.
}

void Spacetime::computeBatchedAder(
    SEISSOL_GPU_PARAM const TimeStepCoefficients& coeffs,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM LocalTmp& tmp,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM recording::ConditionalMaterialTable& materialTable,
    SEISSOL_GPU_PARAM bool updateDisplacement,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(key) == dataTable.end()) {
    return;
  }
  auto& entry = dataTable[key];

  const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
  const auto** localIntegrationPtrs = const_cast<const real**>(
      (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());

  kernel::gpu_derivative derivative = deviceDerivative_;
  kernel::gpu_damageStep step = deviceStep_;
  kernel::gpu_damageTransport assemble = deviceTransport_;
  derivative.numElements = numElements;
  step.numElements = numElements;
  assemble.numElements = numElements;

  const auto maxTmpMem = yateto::getMaxTmpMemRequired(derivative, step, assemble);
  auto tmpMem = runtime.memoryHandle<real>((maxTmpMem * numElements) / sizeof(real));

  // What the recursion transports by, at the cell mean of this step. The
  // geometry it is built from is what the cell keeps, reached the way the
  // step reaches the material: one pointer, an offset further in.
  constexpr auto EpsInitOffset =
      offsetof(LocalIntegrationData, specific) + offsetof(NonLinearLocalData, epsInit);
  constexpr auto ParametersOffset =
      offsetof(LocalIntegrationData, specific) + offsetof(NonLinearLocalData, parameters);
  SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData, starMatrices);
  auto* transportPtrs = (entry.get(inner_keys::Wp::Id::Transport))->getDeviceDataPtr();
  assemble.Q = const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());
  assemble.epsInit = localIntegrationPtrs;
  assemble.extraOffset_epsInit = EpsInitOffset / sizeof(real);
  assemble.materialParameters = localIntegrationPtrs;
  assemble.extraOffset_materialParameters = ParametersOffset / sizeof(real);
  assemble.invariantFloor = invariantFloor();
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    assemble.star(i) = localIntegrationPtrs;
    assemble.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    assemble.transport(i) = transportPtrs;
    assemble.extraOffset_transport(i) = yateto::computeFamilySize<tensor::transport>(1, i);
  }
  assemble.linearAllocator.initialize(tmpMem.get());
  assemble.streamPtr = runtime.stream();
  assemble.execute();

  derivative.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::transport>(); ++i) {
    derivative.transport(i) = const_cast<const real**>(transportPtrs);
    derivative.extraOffset_transport(i) = yateto::computeFamilySize<tensor::transport>(1, i);
  }
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    derivative.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
    derivative.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  derivative.Q =
      const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());
  for (std::size_t der = 0; der < ConvergenceOrder; ++der) {
    derivative.power(der) = coeffs.integral.state[der];
  }
  derivative.linearAllocator.initialize(tmpMem.get());
  derivative.streamPtr = runtime.stream();
  derivative.execute();

  // The rule the step samples with arrives with the step, as it does serially.
  const auto& quadrature = coeffs.quadrature;
  assert(quadrature.nodes.size() == stepNodes());
  for (std::size_t q = 0; q < quadrature.nodes.size(); ++q) {
    for (std::size_t i = 0; i < ConvergenceOrder; ++i) {
      step.evaluate(q, i) = quadrature.coefficients[q].state[i];
    }
    step.weight(q) = quadrature.weights[q];
  }
  step.stepWidth = timeStepWidth;

  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    step.dQ(i) =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr());
    step.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
  }
  // The expansion of what the cell carries beyond its state, behind the
  // expansion of the state in the same buffer -- as it is serially. Without
  // this the step writes the carried columns, the two bounds among them,
  // through a pointer nobody set.
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::transportDer>(); ++i) {
    step.transportDer(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
    step.extraOffset_transportDer(i) = yateto::computeFamilySize<tensor::dQ>() +
                                       yateto::computeFamilySize<tensor::transportDer>(1, i);
  }

  step.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
  step.sourceI = (entry.get(inner_keys::Wp::Id::SourceIntegrals))->getDeviceDataPtr();

  // The initial strain and the material sit in the solver's own part of the
  // cell's local integration data, so they are reached the way an anelastic
  // material's source matrix is: the same pointer, an offset further in.
  constexpr auto EpsInitOffset =
      offsetof(LocalIntegrationData, specific) + offsetof(NonLinearLocalData, epsInit);
  constexpr auto ParametersOffset =
      offsetof(LocalIntegrationData, specific) + offsetof(NonLinearLocalData, parameters);
  static_assert(EpsInitOffset % sizeof(real) == 0 && ParametersOffset % sizeof(real) == 0,
                "The per-cell inputs of the step are not aligned to the real size.");
  step.epsInit = localIntegrationPtrs;
  step.extraOffset_epsInit = EpsInitOffset / sizeof(real);
  step.materialParameters = localIntegrationPtrs;
  step.extraOffset_materialParameters = ParametersOffset / sizeof(real);

  step.invariantFloor = invariantFloor();
  step.linearAllocator.initialize(tmpMem.get());
  step.streamPtr = runtime.stream();
  step.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate Spacetime::metrics() const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::derivative>();
  estimate += PerformanceEstimate::fromKernel<kernel::damageStep>();

  std::uint64_t reals = 0;
  // the state in, the transported tensor and the source integral out
  reals += tensor::Q::size() + tensor::I::size() + tensor::sourceI::size();
  // the star matrices, and the two expansions -- the state's from its first
  // derivative on, because its zeroth is the state and is counted above, and
  // the one the cell carries beyond it, which the step writes in full
  reals += yateto::computeFamilySize<tensor::star>();
  reals += yateto::computeFamilySize<tensor::dQ>(1);
  reals += yateto::computeFamilySize<tensor::transportDer>();

  estimate.bytes = reals * sizeof(real);
  return estimate;
}

void Time::evaluate(const TimeCoefficients& coeffs,
                    const real* timeDerivatives,
                    real timeEvaluated[tensor::I::size()]) {
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeEvaluated)) % Alignment == 0);

  // The columns the transported tensor shares with the state: the Taylor sum
  // the linear solver evaluates, over its own expansion and its own basis.
  kernel::derivativeTaylorExpansion krnl;
  krnl.I = timeEvaluated;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
    krnl.power(i) = coeffs.state[i];
  }
  krnl.execute();

  // And the columns that are not: out of the expansion projected for them,
  // with the coefficients of the basis it was projected into. Reconstructing
  // them rather than rebuilding the stress is the point -- the stress of a
  // cell is a question about that cell's material, and a neighbour has no
  // business answering it.
  kernel::carriedTaylorExpansion carried;
  carried.I = timeEvaluated;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::transportDer>(); ++i) {
    carried.transportDer(i) = timeDerivatives + yateto::computeFamilySize<tensor::dQ>() +
                              yateto::computeFamilySize<tensor::transportDer>(1, i);
    carried.extraPower(i) = coeffs.extra[i];
  }
  carried.execute();
}

void Time::evaluateBatched(SEISSOL_GPU_PARAM const TimeCoefficients& coeffs,
                           SEISSOL_GPU_PARAM const real** timeDerivatives,
                           SEISSOL_GPU_PARAM real** timeIntegratedDofs,
                           SEISSOL_GPU_PARAM std::size_t numElements,
                           SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  assert(timeDerivatives != nullptr);
  assert(timeIntegratedDofs != nullptr);
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

  // Both halves of a reconstruction, as they are serially. The second
  // expansion sits behind the first in the same buffer, so it is the same
  // pointer table an offset further in -- no entry of its own.
  kernel::gpu_derivativeTaylorExpansion krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = timeDerivatives;
    krnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
    krnl.power(i) = coeffs.state[i];
  }
  krnl.streamPtr = runtime.stream();
  krnl.execute();

  kernel::gpu_carriedTaylorExpansion carried;
  carried.numElements = numElements;
  carried.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::transportDer>(); ++i) {
    carried.transportDer(i) = timeDerivatives;
    carried.extraOffset_transportDer(i) = yateto::computeFamilySize<tensor::dQ>() +
                                          yateto::computeFamilySize<tensor::transportDer>(1, i);
    carried.extraPower(i) = coeffs.extra[i];
  }
  carried.streamPtr = runtime.stream();
  carried.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate Time::metrics() const {
  return PerformanceEstimate::fromKernel<kernel::derivativeTaylorExpansion>();
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

void Time::stateToTransport(const real* dofs,
                            const typename seissol::model::MaterialT::Solver::LocalData& local,
                            real* transported) {
  assert(reinterpret_cast<uintptr_t>(dofs) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(transported) % Alignment == 0);

  kernel::stateToTransport krnl;
  krnl.bindGlobals(Pool::host());
  krnl.Q = dofs;
  krnl.I = transported;
  krnl.epsInit = local.epsInit;
  krnl.materialParameters = local.parameters;
  krnl.invariantFloor = invariantFloor();
  krnl.execute();
}

} // namespace seissol::kernels::solver::nonlinearck
