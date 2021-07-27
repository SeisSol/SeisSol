#ifndef EQUATION_DIRICHLET_BOUNDARY_H_
#define EQUATION_DIRICHLET_BOUNDARY_H_

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/typedefs.hpp"

#include "Numerical_aux/Quadrature.h"

namespace {
// Helper functions, needed because C++ doesnt allow partial func. template specialisation  
template<typename MappingKrnl>
void addRotationToProjectKernel(MappingKrnl& projectKernel,
				const CellBoundaryMapping& boundaryMapping) {
  // do nothing
}

//
// GCC warns that the method below is unused. This is not correct.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
template <>
void addRotationToProjectKernel(seissol::kernel::projectToNodalBoundaryRotated& projectKernel,
				const CellBoundaryMapping& boundaryMapping) {
  assert(boundaryMapping.TinvData != nullptr);
  projectKernel.Tinv = boundaryMapping.TinvData;
}
#pragma GCC diagnostic pop

}

namespace seissol {
namespace kernels {

class DirichletBoundary {
 public:

  DirichletBoundary() {
    quadrature::GaussLegendre(quadPoints, quadWeights, CONVERGENCE_ORDER);
  }

  template<typename Func, typename MappingKrnl>
  void evaluate(const real* dofsVolumeInteriorModal,
                int faceIdx,
                const CellBoundaryMapping &boundaryMapping,
                MappingKrnl&& projectKernelPrototype,
                Func&& evaluateBoundaryCondition,
                real* dofsFaceBoundaryNodal) const {
    auto projectKrnl = projectKernelPrototype;
    addRotationToProjectKernel(projectKrnl, boundaryMapping);
    projectKrnl.I = dofsVolumeInteriorModal;
    projectKrnl.INodal = dofsFaceBoundaryNodal;
    projectKrnl.execute(faceIdx);
  
    auto boundaryDofs = init::INodal::view::create(dofsFaceBoundaryNodal);
  
#ifndef MULTIPLE_SIMULATIONS
    static_assert(nodal::tensor::nodes2D::Shape[0] == tensor::INodal::Shape[0],
		  "Need evaluation at all nodes!");
#else
    static_assert(nodal::tensor::nodes2D::Shape[1] == tensor::INodal::Shape[1],
		  "Need evaluation at all nodes!");
#endif

    assert(boundaryMapping.nodes != nullptr);
  
    // Evaluate boundary conditions at precomputed nodes (in global coordinates).
    std::forward<Func>(evaluateBoundaryCondition)(boundaryMapping.nodes, boundaryDofs);
  }

  template<typename Func, typename MappingKrnl>
  void evaluateTimeDependent(const real* dofsVolumeInteriorModal,
			     int faceIdx,
			     const CellBoundaryMapping &boundaryMapping,
			     MappingKrnl&& projectKernelPrototype,
			     Func&& evaluateBoundaryCondition,
			     real* dofsFaceBoundaryNodal,
			     double startTime,
			     double timeStepWidth) const {
    // TODO(Lukas) Implement functions which depend on the interior values...
    auto boundaryDofs = init::INodal::view::create(dofsFaceBoundaryNodal);
  
#ifndef MULTIPLE_SIMULATIONS
    static_assert(nodal::tensor::nodes2D::Shape[0] == tensor::INodal::Shape[0],
		  "Need evaluation at all nodes!");
#else
    static_assert(nodal::tensor::nodes2D::Shape[1] == tensor::INodal::Shape[1],
		  "Need evaluation at all nodes!");
#endif

    assert(boundaryMapping.nodes != nullptr);
  
    // Compute quad points/weights for interval [t, t+dt]
    double timePoints[CONVERGENCE_ORDER];
    double timeWeights[CONVERGENCE_ORDER];
    for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
      timePoints[point] = (timeStepWidth * quadPoints[point] + 2 * startTime + timeStepWidth)/2;
      timeWeights[point] = 0.5 * timeStepWidth * quadWeights[point];
    }
  
    alignas(ALIGNMENT) real dofsFaceBoundaryNodalTmp[tensor::INodal::size()];
    auto boundaryDofsTmp = init::INodal::view::create(dofsFaceBoundaryNodalTmp);
  
    boundaryDofs.setZero();
    boundaryDofsTmp.setZero();
  
    auto updateKernel = kernel::updateINodal{};
    updateKernel.INodal = dofsFaceBoundaryNodal;
    updateKernel.INodalUpdate = dofsFaceBoundaryNodalTmp;
    // Evaluate boundary conditions at precomputed nodes (in global coordinates).
  
    for (int i = 0; i < CONVERGENCE_ORDER; ++i) {
      boundaryDofsTmp.setZero();
      std::forward<Func>(evaluateBoundaryCondition)(boundaryMapping.nodes,
						    timePoints[i],
						    boundaryDofsTmp);
    
      updateKernel.factor = timeWeights[i];
      updateKernel.execute();
    }
  
  }

 private:
  double quadPoints[CONVERGENCE_ORDER];
  double quadWeights[CONVERGENCE_ORDER];
};


void computeAverageDisplacement(double deltaT,
				const real* timeDerivatives,
				const unsigned int derivativesOffsets[CONVERGENCE_ORDER],
				real timeIntegrated[tensor::I::size()]);

} // namespace kernels
} // namespace seissol

#endif
