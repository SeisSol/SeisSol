// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_ANALYTICALBOUNDARY_H_
#define SEISSOL_SRC_KERNELS_ANALYTICALBOUNDARY_H_

#include "Common/Constants.h"
#include "Common/Real.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/ScriptParameters.h"
#include "Numerical/Quadrature.h"
#include "Physics/GhostState.h"
#include "Physics/InitialField.h"
#include "Physics/NonlinearDirichlet.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::kernels {

/**
 * Samples the analytical solution of the scenario at the given nodes and time.
 */
template <typename Cfg>
struct ApplyAnalyticalSolution {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  ApplyAnalyticalSolution(const std::vector<std::unique_ptr<physics::InitialField>>* initConditions,
                          LTS::Ref<Cfg>& data)
      : initConditions_(initConditions), localData_(data) {}

  void operator()(const real* nodes,
                  double time,
                  typename seissol::init::INodal<Cfg>::view::type& boundaryDofs) const {
    assert(initConditions_ != nullptr);

    constexpr auto NodeCount = seissol::tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg>];
    alignas(Alignment) std::array<double, 3> nodesVec[NodeCount];

#pragma omp simd
    for (std::size_t i = 0; i < NodeCount; ++i) {
      nodesVec[i][0] = nodes[i * 3 + 0];
      nodesVec[i][1] = nodes[i * 3 + 1];
      nodesVec[i][2] = nodes[i * 3 + 2];
    }

    // NOTE: not yet tested for multisim setups
    // (only implemented to get the build to work)

    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      auto slicedBoundaryDofs = multisim::simtensor<Cfg>(boundaryDofs, s);
      initConditions_->at(s % initConditions_->size())
          ->evaluate(time,
                     nodesVec,
                     NodeCount,
                     localData_.template get<LTS::Material>(),
                     slicedBoundaryDofs);
    }
  }

  private:
  const std::vector<std::unique_ptr<physics::InitialField>>* initConditions_;
  LTS::Ref<Cfg>& localData_;
};

/**
 * Evaluates the analytical solution at the face nodes and integrates it over the
 * timestep, at one time of the quadrature after the other: for an initial condition that
 * is no script that compiles, which GhostStateBoundary integrates at once.
 */
template <typename Cfg>
class AnalyticalBoundary {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  AnalyticalBoundary() {
    quadrature::GaussLegendre(quadPoints_.data(), quadWeights_.data(), Cfg::ConvergenceOrder);
  }

  template <typename Func>
  void evaluate(const CellBoundaryMapping<Cfg>& boundaryMapping,
                const Func& evaluateBoundaryCondition,
                real* dofsFaceBoundaryNodal,
                double startTime,
                double timeStepWidth) const {
    auto boundaryDofs = init::INodal<Cfg>::view::create(dofsFaceBoundaryNodal);

    static_assert(nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisDim<Cfg>] ==
                      tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg>],
                  "Need evaluation at all nodes!");

    assert(boundaryMapping.nodes != nullptr);

    // Compute quad points/weights for interval [t, t+dt]
    double timePoints[Cfg::ConvergenceOrder];
    double timeWeights[Cfg::ConvergenceOrder];
    for (unsigned point = 0; point < Cfg::ConvergenceOrder; ++point) {
      timePoints[point] = (timeStepWidth * quadPoints_[point] + 2 * startTime + timeStepWidth) / 2;
      timeWeights[point] = 0.5 * timeStepWidth * quadWeights_[point];
    }

    alignas(Alignment) real dofsFaceBoundaryNodalTmp[tensor::INodal<Cfg>::size()];
    auto boundaryDofsTmp = init::INodal<Cfg>::view::create(dofsFaceBoundaryNodalTmp);

    boundaryDofs.setZero();
    boundaryDofsTmp.setZero();

    auto updateKernel = kernel::updateINodal<Cfg>{};
    updateKernel.INodal = dofsFaceBoundaryNodal;
    updateKernel.INodalUpdate = dofsFaceBoundaryNodalTmp;
    // Evaluate boundary conditions at precomputed nodes (in global coordinates).

    for (unsigned i = 0; i < Cfg::ConvergenceOrder; ++i) {
      boundaryDofsTmp.setZero();
      evaluateBoundaryCondition(boundaryMapping.nodes, timePoints[i], boundaryDofsTmp);

      updateKernel.factor = timeWeights[i];
      updateKernel.execute();
    }
  }

  private:
  std::array<double, Cfg::ConvergenceOrder> quadPoints_{};
  std::array<double, Cfg::ConvergenceOrder> quadWeights_{};
};

/**
 * The ghost state of a boundary whose script is integrated over the time step at once
 * (physics::GhostStateIntegral), at the nodes of a face, in global coordinates, into a tensor
 * shaped like INodal, as the nodal flux takes it: the analytical boundary from an initial condition
 * that is a script that compiles, and the nonlinear Dirichlet boundary.
 */
template <typename Cfg>
class GhostStateBoundary {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  static constexpr auto NodeCount = tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg>];
  static constexpr auto ColumnCount = tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  static constexpr auto AlignedColumnCount =
      tensor::INodalAligned<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  static constexpr auto ParameterCount = model::ScriptParameterCount<model::MaterialOf<Cfg>>;

  void bindGlobals(const GlobalData<Cfg>& global) {
    project_.bindGlobals(global);
    projectAligned_.bindGlobals(global);
    rotateToGlobal_.bindGlobals(global);
  }

  /// The analytical boundary at the face `face` of the cell `data`: per simulation, the integral
  /// of the script of its field (`state`) -- or, if the initial condition is no script that
  /// compiles (`state` is null), the fields `initConditions` at one time after the other.
  void analytical(const physics::AnalyticalGhostState* state,
                  const std::vector<std::unique_ptr<physics::InitialField>>* initConditions,
                  LTS::Ref<Cfg>& data,
                  std::size_t face,
                  double start,
                  double width,
                  real* ghostNodal) const {
    const auto& boundaryMapping = data.template get<LTS::BoundaryMapping>()[face];
    if (state == nullptr) {
      assert(initConditions != nullptr);
      timeQuadrature_.evaluate(boundaryMapping,
                               ApplyAnalyticalSolution<Cfg>(initConditions, data),
                               ghostNodal,
                               start,
                               width);
      return;
    }

    const auto points = pointsOf(boundaryMapping);
    const auto parameters = parametersOf(data.template get<LTS::MaterialParameters>());
    auto ghostView = init::INodal<Cfg>::view::create(ghostNodal);
    ghostView.setZero();
    std::array<double, ColumnCount * NodeCount> ghost{};
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      const auto& integral = state->integral(s);
      assert(integral.componentCount() <= ColumnCount);
      integral.evaluate(start,
                        width,
                        state->scriptSimulation(s),
                        points.data(),
                        NodeCount,
                        parameters.data(),
                        nullptr,
                        ghost.data());
      auto ghostOfSimulation = multisim::simtensor<Cfg>(ghostView, s);
      store(ghost.data(), integral.componentCount(), ghostOfSimulation);
    }
  }

  /// The nonlinear Dirichlet boundary at the face `face` of a cell with the parameters
  /// `parameters` (LTS::MaterialParameters): the state inside at the time points of the step, from
  /// `stateAt(tau, state)`, which writes the modal state `tau` into the step to `state` (shaped
  /// like I), taken to the nodes of the face and into the basis of the condition; the integral of
  /// the ghost state taken back into global coordinates.
  template <typename StateAt>
  void nonlinearDirichlet(const physics::NonlinearDirichlet& condition,
                          const StateAt& stateAt,
                          std::size_t face,
                          const CellBoundaryMapping<Cfg>& boundaryMapping,
                          const real* parameters,
                          double start,
                          double width,
                          real* ghostNodal) const {
    const auto& integral = condition.integral();
    const auto parameterValues = parametersOf(parameters);
    const bool faceAligned = condition.faceAligned();
    const auto timePoints = integral.timePointCount();
    const auto innerCount = integral.innerCount();
    assert(integral.componentCount() <= AlignedColumnCount);
    assert(innerCount <= AlignedColumnCount);
    const auto points = pointsOf(boundaryMapping);

    // the state inside per simulation, component k at node p at time point i at
    // [(i * innerCount + k) * NodeCount + p]; none if the script does not read it
    const auto innerSize = timePoints * innerCount * NodeCount;
    thread_local std::vector<double> inner;
    inner.resize(Cfg::NumSimulations * innerSize);
    if (innerCount > 0) {
      thread_local std::vector<double> times;
      times.resize(timePoints);
      integral.timePoints(start, width, times.data());
      for (std::size_t i = 0; i < timePoints; ++i) {
        alignas(Alignment) real state[tensor::I<Cfg>::size()];
        stateAt(times[i] - start, state);
        if (faceAligned) {
          alignas(Alignment) real nodal[tensor::INodalAligned<Cfg>::size()];
          auto project = projectAligned_;
          project.I = state;
          project.Tinv = boundaryMapping.dataTinv;
          project.INodalAligned = nodal;
          project.execute(face);
          load(init::INodalAligned<Cfg>::view::create(nodal), i, innerCount, innerSize, inner);
        } else {
          alignas(Alignment) real nodal[tensor::INodal<Cfg>::size()];
          auto project = project_;
          project.I = state;
          project.INodal = nodal;
          project.execute(face);
          load(init::INodal<Cfg>::view::create(nodal), i, innerCount, innerSize, inner);
        }
      }
    }

    // the integral, in the basis of the condition
    auto ghostView = init::INodal<Cfg>::view::create(ghostNodal);
    ghostView.setZero();
    alignas(Alignment) real aligned[tensor::INodalAligned<Cfg>::size()];
    auto alignedView = init::INodalAligned<Cfg>::view::create(aligned);
    std::array<double, AlignedColumnCount * NodeCount> ghost{};
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      integral.evaluate(start,
                        width,
                        s,
                        points.data(),
                        NodeCount,
                        parameterValues.data(),
                        inner.data() + s * innerSize,
                        ghost.data());
      if (faceAligned) {
        auto alignedOfSimulation = multisim::simtensor<Cfg>(alignedView, s);
        store(ghost.data(), integral.componentCount(), alignedOfSimulation);
      } else {
        auto ghostOfSimulation = multisim::simtensor<Cfg>(ghostView, s);
        store(ghost.data(), integral.componentCount(), ghostOfSimulation);
      }
    }
    if (faceAligned) {
      auto rotate = rotateToGlobal_;
      rotate.INodalAligned = aligned;
      rotate.T = boundaryMapping.dataT;
      rotate.INodal = ghostNodal;
      rotate.execute();
    }
  }

  private:
  /// The parameters of the material of a cell as scripts read them, from their reals.
  static std::array<double, ParameterCount> parametersOf(const real* parameters) {
    std::array<double, ParameterCount> values{};
    for (std::size_t i = 0; i < ParameterCount; ++i) {
      values[i] = parameters[i];
    }
    return values;
  }

  /// The nodes of the face, as a script reads them.
  static std::array<std::array<double, 3>, NodeCount>
      pointsOf(const CellBoundaryMapping<Cfg>& boundaryMapping) {
    assert(boundaryMapping.nodes != nullptr);
    std::array<std::array<double, 3>, NodeCount> points{};
    for (std::size_t p = 0; p < NodeCount; ++p) {
      for (std::size_t d = 0; d < 3; ++d) {
        points[p][d] = boundaryMapping.nodes[p * 3 + d];
      }
    }
    return points;
  }

  /// Component j at node p of `values`, at values[j * NodeCount + p], into column j of `target`.
  template <typename ViewT>
  static void store(const double* values, std::size_t components, ViewT& target) {
    for (std::size_t j = 0; j < components; ++j) {
      for (std::size_t p = 0; p < NodeCount; ++p) {
        target(p, j) = static_cast<real>(values[j * NodeCount + p]);
      }
    }
  }

  /// The first `innerCount` columns of the state `nodal` at the nodes, of all simulations, as the
  /// time point `point` of `inner`.
  template <typename ViewT>
  static void load(ViewT nodal,
                   std::size_t point,
                   std::size_t innerCount,
                   std::size_t innerSize,
                   std::vector<double>& inner) {
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      auto nodalOfSimulation = multisim::simtensor<Cfg>(nodal, s);
      double* target = inner.data() + s * innerSize + point * innerCount * NodeCount;
      for (std::size_t k = 0; k < innerCount; ++k) {
        for (std::size_t p = 0; p < NodeCount; ++p) {
          target[k * NodeCount + p] = nodalOfSimulation(p, k);
        }
      }
    }
  }

  AnalyticalBoundary<Cfg> timeQuadrature_;
  kernel::projectToFaceNodes<Cfg> project_;
  kernel::projectToFaceNodesAligned<Cfg> projectAligned_;
  kernel::rotateNodalToGlobal<Cfg> rotateToGlobal_;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_ANALYTICALBOUNDARY_H_
