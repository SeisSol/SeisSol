// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_ANALYTICALBOUNDARY_H_
#define SEISSOL_SRC_KERNELS_ANALYTICALBOUNDARY_H_

#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Numerical/Quadrature.h"
#include "Physics/InitialField.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cassert>
#include <memory>
#include <vector>

namespace seissol::kernels {

/**
 * Samples the analytical solution of the scenario at the given nodes and time.
 */
struct ApplyAnalyticalSolution {
  ApplyAnalyticalSolution(const std::vector<std::unique_ptr<physics::InitialField>>* initConditions,
                          LTS::Ref& data)
      : initConditions_(initConditions), localData_(data) {}

  void operator()(const real* nodes,
                  double time,
                  seissol::init::INodal::view::type& boundaryDofs) const {
    assert(initConditions_ != nullptr);

    constexpr auto NodeCount = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension];
    alignas(Alignment) std::array<double, 3> nodesVec[NodeCount];

#pragma omp simd
    for (std::size_t i = 0; i < NodeCount; ++i) {
      nodesVec[i][0] = nodes[i * 3 + 0];
      nodesVec[i][1] = nodes[i * 3 + 1];
      nodesVec[i][2] = nodes[i * 3 + 2];
    }

    // NOTE: not yet tested for multisim setups
    // (only implemented to get the build to work)

    for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
      auto slicedBoundaryDofs = multisim::simtensor(boundaryDofs, s);
      initConditions_->at(s % initConditions_->size())
          ->evaluate(
              time, nodesVec, NodeCount, localData_.get<LTS::Material>(), slicedBoundaryDofs);
    }
  }

  private:
  const std::vector<std::unique_ptr<physics::InitialField>>* initConditions_;
  LTS::Ref& localData_;
};

/**
 * Evaluates the analytical solution at the face nodes and integrates it over the
 * timestep. The condition is a function of position and time alone; one that also
 * depended on the interior state would need the time-integrated DOFs passed in.
 */
class AnalyticalBoundary {
  public:
  AnalyticalBoundary() { quadrature::GaussLegendre(quadPoints_, quadWeights_, ConvergenceOrder); }

  template <typename Func>
  void evaluate(const CellBoundaryMapping& boundaryMapping,
                const Func& evaluateBoundaryCondition,
                real* dofsFaceBoundaryNodal,
                double startTime,
                double timeStepWidth) const {
    auto boundaryDofs = init::INodal::view::create(dofsFaceBoundaryNodal);

    static_assert(nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension] ==
                      tensor::INodal::Shape[multisim::BasisFunctionDimension],
                  "Need evaluation at all nodes!");

    assert(boundaryMapping.nodes != nullptr);

    // Compute quad points/weights for interval [t, t+dt]
    double timePoints[ConvergenceOrder];
    double timeWeights[ConvergenceOrder];
    for (unsigned point = 0; point < ConvergenceOrder; ++point) {
      timePoints[point] = (timeStepWidth * quadPoints_[point] + 2 * startTime + timeStepWidth) / 2;
      timeWeights[point] = 0.5 * timeStepWidth * quadWeights_[point];
    }

    alignas(Alignment) real dofsFaceBoundaryNodalTmp[tensor::INodal::size()];
    auto boundaryDofsTmp = init::INodal::view::create(dofsFaceBoundaryNodalTmp);

    boundaryDofs.setZero();
    boundaryDofsTmp.setZero();

    auto updateKernel = kernel::updateINodal{};
    updateKernel.INodal = dofsFaceBoundaryNodal;
    updateKernel.INodalUpdate = dofsFaceBoundaryNodalTmp;
    // Evaluate boundary conditions at precomputed nodes (in global coordinates).

    for (unsigned i = 0; i < ConvergenceOrder; ++i) {
      boundaryDofsTmp.setZero();
      evaluateBoundaryCondition(boundaryMapping.nodes, timePoints[i], boundaryDofsTmp);

      updateKernel.factor = timeWeights[i];
      updateKernel.execute();
    }
  }

  private:
  double quadPoints_[ConvergenceOrder]{};
  double quadWeights_[ConvergenceOrder]{};
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_ANALYTICALBOUNDARY_H_
