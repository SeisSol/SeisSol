#ifndef EQUATION_DIRICHLET_BOUNDARY_H_
#define EQUATION_DIRICHLET_BOUNDARY_H_

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/typedefs.hpp"

#include "Numerical_aux/Quadrature.h"

#ifdef ACL_DEVICE
#include "yateto.h"
#include "device.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.hpp"
#include "Equations/elastic/Kernels/DeviceAux/KernelsAux.h"
#endif

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
  
    static_assert(nodal::tensor::nodes2D::Shape[0] == tensor::INodal::Shape[0],
		  "Need evaluation at all nodes!");

    assert(boundaryMapping.nodes != nullptr);
  
    // Evaluate boundary conditions at precomputed nodes (in global coordinates).
    std::forward<Func>(evaluateBoundaryCondition)(boundaryMapping.nodes, boundaryDofs);
  }

#ifdef ACL_DEVICE
  template<typename Func, typename MappingKrnl, typename InverseMappingKrnl>
  void evaluateOnDevice(int faceIdx,
                        ConditionalKey& key,
                        MappingKrnl& projectKernelPrototype,
                        InverseMappingKrnl& nodalLfKrnlPrototype,
                        local_flux::aux::DirichletBoundaryAux<Func>& boundaryCondition,
                        ConditionalPointersToRealsTable &dataTable,
                        device::DeviceInstance& device) const {

    const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Dofs)->getSize()};

    size_t memCounter{0};
    auto* dofsFaceBoundaryNodalData = reinterpret_cast<real*>(device.api->getStackMemory(tensor::INodal::size() * numElements * sizeof(real)));
    auto** dofsFaceBoundaryNodalPtrs = reinterpret_cast<real**>(device.api->getStackMemory(numElements * sizeof(real*)));
    memCounter += 2;

    auto* deviceStream = device.api->getDefaultStream();
    device.algorithms.incrementalAdd(
      dofsFaceBoundaryNodalPtrs,
      dofsFaceBoundaryNodalData,
      tensor::INodal::size(),
      numElements,
      deviceStream
    );

    const auto auxTmpMemSize = yateto::getMaxTmpMemRequired(nodalLfKrnlPrototype, projectKernelPrototype);
    auto* auxTmpMem = reinterpret_cast<real*>(device.api->getStackMemory(auxTmpMemSize * numElements));
    memCounter += 1;

    auto** TinvData = dataTable[key].get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();
    auto** idofsPtrs = dataTable[key].get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr();

    auto projectKrnl = projectKernelPrototype;
    projectKrnl.numElements = numElements;
    projectKrnl.Tinv = const_cast<const real **>(TinvData);
    projectKrnl.I = const_cast<const real **>(idofsPtrs);
    projectKrnl.INodal = dofsFaceBoundaryNodalPtrs;
    projectKrnl.linearAllocator.initialize(auxTmpMem);
    projectKrnl.streamPtr = deviceStream;
    projectKrnl.execute(faceIdx);


    boundaryCondition.evaluate(dofsFaceBoundaryNodalPtrs, numElements, deviceStream);


    auto** dofsPtrs = dataTable[key].get(inner_keys::Wp::Id::Dofs)->getDeviceDataPtr();
    auto** nAmNm1 = dataTable[key].get(inner_keys::Wp::Id::AminusT)->getDeviceDataPtr();

    auto nodalLfKrnl = nodalLfKrnlPrototype;
    nodalLfKrnl.numElements = numElements;
    nodalLfKrnl.Q = dofsPtrs;
    nodalLfKrnl.INodal = const_cast<const real **>(dofsFaceBoundaryNodalPtrs);
    nodalLfKrnl.AminusT = const_cast<const real **>(nAmNm1);
    nodalLfKrnl.linearAllocator.initialize(auxTmpMem);
    nodalLfKrnl.streamPtr = deviceStream;
    nodalLfKrnl.execute(faceIdx);

    for (size_t i = 0; i < memCounter; ++i)
      device.api->popStackMemory();
  }
#endif

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
  
    static_assert(nodal::tensor::nodes2D::Shape[0] == tensor::INodal::Shape[0],
		  "Need evaluation at all nodes!");

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
