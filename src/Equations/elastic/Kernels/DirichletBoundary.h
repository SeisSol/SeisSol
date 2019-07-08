#ifndef EQUATION_DIRICHLET_BOUNDARY_H_
#define EQUATION_DIRICHLET_BOUNDARY_H_

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"
#include "Initializer/typedefs.hpp"

namespace {
// Helper functions, needed because C++ doesnt allow partial func. template specialisation  
template<typename MappingKrnl>
void addRotationToProjectKernel(MappingKrnl& projectKernelPrototype,
				const CellBoundaryMapping& boundaryMapping) {
  // do nothing
}
 
template <>
void addRotationToProjectKernel(seissol::kernel::projectToNodalBoundaryRotated& projectKernelPrototype,
				const CellBoundaryMapping& boundaryMapping) {
  projectKernelPrototype.T = boundaryMapping.TData;
}

}

namespace seissol {

template<typename Func, typename MappingKrnl>
void computeDirichletBoundary(const real* dofsVolumeInteriorModal,
				     int faceIdx,
				     const CellBoundaryMapping &boundaryMapping,
				     MappingKrnl&& projectKernelPrototype,
				     Func&& evaluateBoundaryCondition,
				     real* dofsFaceBoundaryNodal) {
  // TODO(Lukas) Implement time-dependent functions
  
  auto projectKrnl = projectKernelPrototype;
  addRotationToProjectKernel(projectKernelPrototype, boundaryMapping);
  projectKrnl.I = dofsVolumeInteriorModal;
  projectKrnl.INodal = dofsFaceBoundaryNodal;
  projectKrnl.execute(faceIdx);
  
  auto boundaryDofs = init::INodal::view::create(dofsFaceBoundaryNodal);
  
  // Evaluate boundary conditions at precomputed nodes (in global coordinates).
  const auto &nodesData = boundaryMapping.nodes;
  auto nodes = init::nodes2D::view::create(nodesData);
  
  static_assert(tensor::nodes2D::Shape[0] == tensor::INodal::Shape[0],
		"Need evaluation at all nodes!");
  
  // Set values for all nodes
  std::forward<Func>(evaluateBoundaryCondition)(nodes, boundaryDofs);
}

void computeAverageDisplacement(double deltaT,
					 const real* timeDerivatives,
					 const unsigned int derivativesOffsets[CONVERGENCE_ORDER],
					 real timeIntegrated[tensor::I::size()]);

} // namespace seissol

#endif
