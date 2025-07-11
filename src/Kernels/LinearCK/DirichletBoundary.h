// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DIRICHLETBOUNDARY_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DIRICHLETBOUNDARY_H_

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/Typedefs.h"

#include "Common/Offset.h"

#include "Numerical/Quadrature.h"
#include <Parallel/Runtime/Stream.h>

#include "Solver/MultipleSimulations.h"

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/LinearCK/DeviceAux/KernelsAux.h"
#include "device.h"
#include "yateto.h"
#endif

namespace {
// Helper functions, needed because C++ doesnt allow partial func. template specialisation
template <typename MappingKrnl>
void addRotationToProjectKernel(MappingKrnl& projectKernel,
                                const seissol::CellBoundaryMapping& boundaryMapping) {
  // do nothing
}

//
// GCC warns that the method below is unused. This is not correct.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
template <>
void addRotationToProjectKernel(seissol::kernel::projectToNodalBoundaryRotated& projectKernel,
                                const seissol::CellBoundaryMapping& boundaryMapping) {
  assert(boundaryMapping.dataTinv != nullptr);
  projectKernel.Tinv = boundaryMapping.dataTinv;
}
#pragma GCC diagnostic pop

} // namespace

namespace seissol::kernels {

class DirichletBoundary {
  public:
  DirichletBoundary() { quadrature::GaussLegendre(quadPoints, quadWeights, ConvergenceOrder); }

  template <typename Func, typename MappingKrnl>
  void evaluate(const real* dofsVolumeInteriorModal,
                int faceIdx,
                const CellBoundaryMapping& boundaryMapping,
                MappingKrnl&& projectKernelPrototype,
                Func&& evaluateBoundaryCondition,
                real* dofsFaceBoundaryNodal) const {
    auto projectKrnl = std::forward<MappingKrnl>(projectKernelPrototype);
    addRotationToProjectKernel(projectKrnl, boundaryMapping);
    projectKrnl.I = dofsVolumeInteriorModal;
    projectKrnl.INodal = dofsFaceBoundaryNodal;
    projectKrnl.execute(faceIdx);

    auto boundaryDofs = init::INodal::view::create(dofsFaceBoundaryNodal);

    static_assert(nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension] ==
                      tensor::INodal::Shape[multisim::BasisFunctionDimension],
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
                        device::DeviceInstance& device,
                        seissol::parallel::runtime::StreamRuntime& runtime) const {

    const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Dofs)->getSize()};

    auto** dofsFaceBoundaryNodalPtrs =
        dataTable[key].get(inner_keys::Wp::Id::DofsFaceBoundaryNodal)->getDeviceDataPtr();

    auto* deviceStream = runtime.stream();

    const auto auxTmpMemSize =
        yateto::getMaxTmpMemRequired(nodalLfKrnlPrototype, projectKernelPrototype);
    auto auxTmpMem = runtime.memoryHandle<real>((auxTmpMemSize * numElements) / sizeof(real));

    auto** dataTinv = dataTable[key].get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();
    auto** idofsPtrs = dataTable[key].get(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr();

    auto projectKrnl = projectKernelPrototype;
    projectKrnl.numElements = numElements;
    projectKrnl.Tinv = const_cast<const real**>(dataTinv);
    projectKrnl.I = const_cast<const real**>(idofsPtrs);
    projectKrnl.INodal = dofsFaceBoundaryNodalPtrs;
    projectKrnl.linearAllocator.initialize(auxTmpMem.get());
    projectKrnl.streamPtr = deviceStream;
    projectKrnl.execute(faceIdx);

    boundaryCondition.evaluate(dofsFaceBoundaryNodalPtrs, numElements, deviceStream);

    auto** dofsPtrs = dataTable[key].get(inner_keys::Wp::Id::Dofs)->getDeviceDataPtr();

    auto nodalLfKrnl = nodalLfKrnlPrototype;
    nodalLfKrnl.numElements = numElements;
    nodalLfKrnl.Q = dofsPtrs;
    nodalLfKrnl.INodal = const_cast<const real**>(dofsFaceBoundaryNodalPtrs);
    nodalLfKrnl.AminusT = const_cast<const real**>(
        dataTable[key].get(inner_keys::Wp::Id::NeighborIntegrationData)->getDeviceDataPtr());
    nodalLfKrnl.extraOffset_AminusT =
        SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, faceIdx);
    nodalLfKrnl.linearAllocator.initialize(auxTmpMem.get());
    nodalLfKrnl.streamPtr = deviceStream;
    nodalLfKrnl.execute(faceIdx);
  }
#endif

  template <typename Func, typename MappingKrnl>
  void evaluateTimeDependent(const real* dofsVolumeInteriorModal,
                             int faceIdx,
                             const CellBoundaryMapping& boundaryMapping,
                             const MappingKrnl& projectKernelPrototype,
                             Func&& evaluateBoundaryCondition,
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
      std::forward<Func>(evaluateBoundaryCondition)(
          boundaryMapping.nodes, timePoints[i], boundaryDofsTmp);

      updateKernel.factor = timeWeights[i];
      updateKernel.execute();
    }
  }

  private:
  double quadPoints[ConvergenceOrder];
  double quadWeights[ConvergenceOrder];
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DIRICHLETBOUNDARY_H_
