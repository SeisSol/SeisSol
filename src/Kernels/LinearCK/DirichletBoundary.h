// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DIRICHLETBOUNDARY_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DIRICHLETBOUNDARY_H_

#include "Common/Offset.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Numerical/Quadrature.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/LinearCK/DeviceAux/KernelsAux.h"

#include <Device/device.h>
#include <yateto.h>
#endif

namespace seissol::kernels {

class DirichletBoundary {
  public:
  DirichletBoundary() { quadrature::GaussLegendre(quadPoints, quadWeights, ConvergenceOrder); }

  template <typename Func, typename MappingKrnl>
  void evaluateTimeDependent(const real* /*dofsVolumeInteriorModal*/,
                             int /*faceIdx*/,
                             const CellBoundaryMapping& boundaryMapping,
                             const MappingKrnl& /*projectKernelPrototype*/,
                             const Func& evaluateBoundaryCondition,
                             real* dofsFaceBoundaryNodal,
                             double startTime,
                             double timeStepWidth) const {
    // TODO(Lukas) Implement functions which depend on the interior values...
    auto boundaryDofs = init::INodal::view::create(dofsFaceBoundaryNodal);

    static_assert(nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension] ==
                      tensor::INodal::Shape[multisim::BasisFunctionDimension],
                  "Need evaluation at all nodes!");

    assert(boundaryMapping.nodes != nullptr);

    // Compute quad points/weights for interval [t, t+dt]
    double timePoints[ConvergenceOrder];
    double timeWeights[ConvergenceOrder];
    for (unsigned point = 0; point < ConvergenceOrder; ++point) {
      timePoints[point] = (timeStepWidth * quadPoints[point] + 2 * startTime + timeStepWidth) / 2;
      timeWeights[point] = 0.5 * timeStepWidth * quadWeights[point];
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
  double quadPoints[ConvergenceOrder]{};
  double quadWeights[ConvergenceOrder]{};
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DIRICHLETBOUNDARY_H_
