#ifndef EQUATION_DIRICHLET_BOUNDARY_H_
#define EQUATION_DIRICHLET_BOUNDARY_H_

#include "Initializer/typedefs.hpp"

#include "Numerical_aux/Quadrature.h"
#include <Common/configtensor.hpp>

#ifdef ACL_DEVICE
#include "yateto.h"
#include "device.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.hpp"
#include "Equations/elastic/Kernels/DeviceAux/KernelsAux.h"
#endif

namespace {
// Helper functions, needed because C++ doesnt allow partial func. template specialisation
template <typename MappingKrnl>
void addRotationToProjectKernel(MappingKrnl& projectKernel,
                                const CellBoundaryMapping& boundaryMapping) {
  // do nothing
}

//
// GCC warns that the method below is unused. This is not correct.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
template <>
void addRotationToProjectKernel(
    seissol::Yateto<Config>::Kernel::projectToNodalBoundaryRotated& projectKernel,
    const CellBoundaryMapping& boundaryMapping) {
  assert(boundaryMapping.TinvData != nullptr);
  projectKernel.Tinv = boundaryMapping.TinvData;
}
#pragma GCC diagnostic pop

} // namespace

namespace seissol::kernels {

template <typename Config>
class DirichletBoundary {
  public:
  using RealT = typename Config::RealT;

  DirichletBoundary() {
    quadrature::GaussLegendre(quadPoints, quadWeights, Config::ConvergenceOrder);
  }

  template <typename Func, typename MappingKrnl>
  void evaluate(const RealT* dofsVolumeInteriorModal,
                int faceIdx,
                const CellBoundaryMapping& boundaryMapping,
                MappingKrnl&& projectKernelPrototype,
                Func&& evaluateBoundaryCondition,
                RealT* dofsFaceBoundaryNodal) const {
    auto projectKrnl = projectKernelPrototype;
    addRotationToProjectKernel(projectKrnl, boundaryMapping);
    projectKrnl.I = dofsVolumeInteriorModal;
    projectKrnl.INodal = dofsFaceBoundaryNodal;
    projectKrnl.execute(faceIdx);

    auto boundaryDofs = Yateto<Config>::Init::INodal::view::create(dofsFaceBoundaryNodal);

    static_assert(Yateto<Config>::Tensor::nodal::nodes2D::Shape[0] ==
                      Yateto<Config>::Tensor::INodal::Shape[0],
                  "Need evaluation at all nodes!");

    assert(boundaryMapping.nodes != nullptr);

    // Evaluate boundary conditions at precomputed nodes (in global coordinates).
    std::forward<Func>(evaluateBoundaryCondition)(boundaryMapping.nodes, boundaryDofs);
  }

#ifdef ACL_DEVICE
  template <typename Func, typename MappingKrnl, typename InverseMappingKrnl>
  void evaluateOnDevice(int faceIdx,
                        ConditionalKey& key,
                        MappingKrnl& projectKernelPrototype,
                        InverseMappingKrnl& nodalLfKrnlPrototype,
                        local_flux::aux::DirichletBoundaryAux<Func>& boundaryCondition,
                        ConditionalPointersToRealsTable& dataTable,
                        device::DeviceInstance& device) const {

    const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Dofs)->getSize()};

    size_t memCounter{0};
    auto* dofsFaceBoundaryNodalData = reinterpret_cast<real*>(device.api->getStackMemory(
        Yateto<Config>::Tensor::INodal::size() * numElements * sizeof(real)));
    auto** dofsFaceBoundaryNodalPtrs =
        reinterpret_cast<real**>(device.api->getStackMemory(numElements * sizeof(real*)));
    memCounter += 2;

    auto* deviceStream = device.api->getDefaultStream();
    device.algorithms.incrementalAdd(dofsFaceBoundaryNodalPtrs,
                                     dofsFaceBoundaryNodalData,
                                     Yateto<Config>::Tensor::INodal::size(),
                                     numElements,
                                     deviceStream);

    const auto auxTmpMemSize =
        yateto::getMaxTmpMemRequired(nodalLfKrnlPrototype, projectKernelPrototype);
    auto* auxTmpMem =
        reinterpret_cast<real*>(device.api->getStackMemory(auxTmpMemSize * numElements));
    memCounter += 1;

    auto** TinvData = dataTable[key].get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();
    auto** idofsPtrs = dataTable[key].get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr();

    auto projectKrnl = projectKernelPrototype;
    projectKrnl.numElements = numElements;
    projectKrnl.Tinv = const_cast<const real**>(TinvData);
    projectKrnl.I = const_cast<const real**>(idofsPtrs);
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
    nodalLfKrnl.INodal = const_cast<const real**>(dofsFaceBoundaryNodalPtrs);
    nodalLfKrnl.AminusT = const_cast<const real**>(nAmNm1);
    nodalLfKrnl.linearAllocator.initialize(auxTmpMem);
    nodalLfKrnl.streamPtr = deviceStream;
    nodalLfKrnl.execute(faceIdx);

    for (size_t i = 0; i < memCounter; ++i)
      device.api->popStackMemory();
  }
#endif

  template <typename Func, typename MappingKrnl>
  void evaluateTimeDependent(const RealT* dofsVolumeInteriorModal,
                             int faceIdx,
                             const CellBoundaryMapping& boundaryMapping,
                             MappingKrnl&& projectKernelPrototype,
                             Func&& evaluateBoundaryCondition,
                             RealT* dofsFaceBoundaryNodal,
                             double startTime,
                             double timeStepWidth) const {
    // TODO(Lukas) Implement functions which depend on the interior values...
    auto boundaryDofs = Yateto<Config>::Init::INodal::view::create(dofsFaceBoundaryNodal);

    static_assert(Yateto<Config>::Tensor::nodal::nodes2D::Shape[0] ==
                      Yateto<Config>::Tensor::INodal::Shape[0],
                  "Need evaluation at all nodes!");

    assert(boundaryMapping.nodes != nullptr);

    // Compute quad points/weights for interval [t, t+dt]
    double timePoints[Config::ConvergenceOrder];
    double timeWeights[Config::ConvergenceOrder];
    for (unsigned point = 0; point < Config::ConvergenceOrder; ++point) {
      timePoints[point] = (timeStepWidth * quadPoints[point] + 2 * startTime + timeStepWidth) / 2;
      timeWeights[point] = 0.5 * timeStepWidth * quadWeights[point];
    }

    alignas(Alignment) RealT dofsFaceBoundaryNodalTmp[Yateto<Config>::Tensor::INodal::size()];
    auto boundaryDofsTmp = Yateto<Config>::Init::INodal::view::create(dofsFaceBoundaryNodalTmp);

    boundaryDofs.setZero();
    boundaryDofsTmp.setZero();

    auto updateKernel = typename Yateto<Config>::Kernel::updateINodal{};
    updateKernel.INodal = dofsFaceBoundaryNodal;
    updateKernel.INodalUpdate = dofsFaceBoundaryNodalTmp;
    // Evaluate boundary conditions at precomputed nodes (in global coordinates).

    for (int i = 0; i < Config::ConvergenceOrder; ++i) {
      boundaryDofsTmp.setZero();
      std::forward<Func>(evaluateBoundaryCondition)(
          boundaryMapping.nodes, timePoints[i], boundaryDofsTmp);

      updateKernel.factor = timeWeights[i];
      updateKernel.execute();
    }
  }

  private:
  double quadPoints[Config::ConvergenceOrder];
  double quadWeights[Config::ConvergenceOrder];
};

template <typename Config>
void computeAverageDisplacement(
    double deltaT,
    const typename Config::RealT* timeDerivatives,
    const unsigned int derivativesOffsets[Config::ConvergenceOrder],
    typename Config::RealT timeIntegrated[ConfigConstants<Config>::TensorSizeI]) {
  // TODO(Lukas) Only compute integral for displacement, not for all vars.
  assert(reinterpret_cast<uintptr_t>(timeDerivatives) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(timeIntegrated) % Alignment == 0);
  assert(deltaT > 0);

  typename Yateto<Config>::Kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeIntegrated;
  for (size_t i = 0; i < Config::ConvergenceOrder; ++i) {
    intKrnl.dQ(i) = timeDerivatives + derivativesOffsets[i];
  }

  typename Config::RealT factorial = 2.0;
  double power = deltaT * deltaT;

  for (int der = 0; der < Config::ConvergenceOrder; ++der) {
    intKrnl.power = power / factorial;
    intKrnl.execute(der);

    factorial *= der + 2.0;
    power *= deltaT;
  }
}

} // namespace seissol::kernels

#endif
