// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Kernels/Local.h"

#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/LocalBase.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Parallel/Runtime/Stream.h>
#include <Physics/InitialField.h>
#include <cstddef>
#include <generated_code/init.h>
#include <generated_code/kernel.h>
#include <tensor.h>
#include <vector>
#include <yateto.h>

#include <array>
#include <cassert>
#include <stdint.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop

#include "Kernels/Common.h"

#include "utils/logger.h"

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)
namespace seissol::kernels {

void LocalBase::checkGlobalData(const GlobalData* global, size_t alignment) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert((reinterpret_cast<uintptr_t>(global->stiffnessMatrices(stiffness))) % alignment == 0);
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert((reinterpret_cast<uintptr_t>(global->localChangeOfBasisMatricesTransposed(flux))) %
               alignment ==
           0);
    assert((reinterpret_cast<uintptr_t>(global->changeOfBasisMatrices(flux))) % alignment == 0);
  }
#endif
}

void Local::setHostGlobalData(const GlobalData* global) {
  checkGlobalData(global, Alignment);
  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;

  m_nodalLfKrnlPrototype.project2nFaceTo3m = global->project2nFaceTo3m;

  m_projectKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
  m_projectRotatedKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
}

void Local::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);

  deviceVolumeKernelPrototype.kDivM = global.onDevice->stiffnessMatrices;
#ifdef USE_PREMULTIPLY_FLUX
  deviceLocalFluxKernelPrototype.plusFluxMatrices = global.onDevice->plusFluxMatrices;
#else
  deviceLocalFluxKernelPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceLocalFluxKernelPrototype.fMrT = global.onDevice->localChangeOfBasisMatricesTransposed;
#endif
  deviceNodalLfKrnlPrototype.project2nFaceTo3m = global.onDevice->project2nFaceTo3m;
  deviceProjectRotatedKrnlPrototype.V3mTo2nFace = global.onDevice->V3mTo2nFace;
#endif
}

template <typename LocalDataType>
struct ApplyAnalyticalSolution {
  ApplyAnalyticalSolution(seissol::physics::InitialField* initCondition, LocalDataType& data)
      : initCondition(initCondition), localData(data) {}

  void operator()(const real* nodes, double time, seissol::init::INodal::view::type& boundaryDofs) {

    auto nodesVec = std::vector<std::array<double, 3>>{};
    int offset = 0;
    for (unsigned int i = 0; i < seissol::tensor::INodal::Shape[0]; ++i) {
      auto curNode = std::array<double, 3>{};
      curNode[0] = nodes[offset++];
      curNode[1] = nodes[offset++];
      curNode[2] = nodes[offset++];
      nodesVec.push_back(curNode);
    }

    assert(initCondition != nullptr);
    initCondition->evaluate(time, nodesVec, localData.material(), boundaryDofs);
  }

  private:
  seissol::physics::InitialField* initCondition{};
  LocalDataType& localData;
};

void Local::computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                            LocalData& data,
                            LocalTmp& tmp,
                            // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                            const CellMaterialData* materialData,
                            const CellBoundaryMapping (*cellBoundaryMapping)[4],
                            double time,
                            double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(timeIntegratedDegreesOfFreedom) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(data.dofs()) % Alignment == 0);

  kernel::volume volKrnl = m_volumeKernelPrototype;
  volKrnl.Q = data.dofs();
  volKrnl.I = timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration().starMatrices[i];
  }

  // Optional source term
  set_ET(volKrnl, get_ptr_sourceMatrix(data.localIntegration().specific));

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.dofs();
  lfKrnl.I = timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs() + tensor::Q::size();

  volKrnl.execute();

  for (int face = 0; face < 4; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.cellInformation().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.localIntegration().nApNm1[face];
      lfKrnl.execute(face);
    }

    alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal::size()];
    auto nodalLfKrnl = m_nodalLfKrnlPrototype;
    nodalLfKrnl.Q = data.dofs();
    nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
    nodalLfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I::size();
    nodalLfKrnl._prefetch.Q = data.dofs() + tensor::Q::size();
    nodalLfKrnl.AminusT = data.neighboringIntegration().nAmNm1[face];

    // Include some boundary conditions here.
    switch (data.cellInformation().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      assert(cellBoundaryMapping != nullptr);
      assert(materialData != nullptr);
      auto* displ = tmp.nodalAvgDisplacements[face].data();
      auto displacement = init::averageNormalDisplacement::view::create(displ);
      // lambdas can't catch gravitationalAcceleration directly, so have to make a copy here.
      const auto localG = gravitationalAcceleration;
      auto applyFreeSurfaceBc =
          [&displacement, &materialData, &localG](const real*, // nodes are unused
                                                  init::INodal::view::type& boundaryDofs) {
            for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
              const double rho = materialData->local.rho;
              assert(localG > 0);
              const double pressureAtBnd = -1 * rho * localG * displacement(i);

              boundaryDofs(i, 0) = 2 * pressureAtBnd - boundaryDofs(i, 0);
              boundaryDofs(i, 1) = 2 * pressureAtBnd - boundaryDofs(i, 1);
              boundaryDofs(i, 2) = 2 * pressureAtBnd - boundaryDofs(i, 2);
            }
          };

      dirichletBoundary.evaluate(timeIntegratedDegreesOfFreedom,
                                 face,
                                 (*cellBoundaryMapping)[face],
                                 m_projectRotatedKrnlPrototype,
                                 applyFreeSurfaceBc,
                                 dofsFaceBoundaryNodal);

      nodalLfKrnl.execute(face);
      break;
    }
    case FaceType::Dirichlet: {
      assert(cellBoundaryMapping != nullptr);
      auto* easiBoundaryMap = (*cellBoundaryMapping)[face].easiBoundaryMap;
      auto* easiBoundaryConstant = (*cellBoundaryMapping)[face].easiBoundaryConstant;
      assert(easiBoundaryConstant != nullptr);
      assert(easiBoundaryMap != nullptr);
      auto applyEasiBoundary = [easiBoundaryMap, easiBoundaryConstant](
                                   const real* nodes, init::INodal::view::type& boundaryDofs) {
        seissol::kernel::createEasiBoundaryGhostCells easiBoundaryKernel;
        easiBoundaryKernel.easiBoundaryMap = easiBoundaryMap;
        easiBoundaryKernel.easiBoundaryConstant = easiBoundaryConstant;
        easiBoundaryKernel.easiIdentMap = init::easiIdentMap::Values;
        easiBoundaryKernel.INodal = boundaryDofs.data();
        easiBoundaryKernel.execute();
      };

      // Compute boundary in [n, t_1, t_2] basis
      dirichletBoundary.evaluate(timeIntegratedDegreesOfFreedom,
                                 face,
                                 (*cellBoundaryMapping)[face],
                                 m_projectRotatedKrnlPrototype,
                                 applyEasiBoundary,
                                 dofsFaceBoundaryNodal);

      // We do not need to rotate the boundary data back to the [x,y,z] basis
      // as we set the Tinv matrix to the identity matrix in the flux solver
      // See init. in CellLocalMatrices.initializeCellLocalMatrices!

      nodalLfKrnl.execute(face);
      break;
    }
    case FaceType::Analytical: {
      assert(cellBoundaryMapping != nullptr);
      assert(initConds != nullptr);
      assert(initConds->size() == 1);
      ApplyAnalyticalSolution applyAnalyticalSolution(this->getInitCond(0), data);

      dirichletBoundary.evaluateTimeDependent(timeIntegratedDegreesOfFreedom,
                                              face,
                                              (*cellBoundaryMapping)[face],
                                              m_projectKrnlPrototype,
                                              applyAnalyticalSolution,
                                              dofsFaceBoundaryNodal,
                                              time,
                                              timeStepWidth);
      nodalLfKrnl.execute(face);
      break;
    }
    default:
      // No boundary condition.
      break;
    }
  }
}

void Local::computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                   ConditionalMaterialTable& materialTable,
                                   ConditionalIndicesTable& indicesTable,
                                   kernels::LocalData::Loader& loader,
                                   LocalTmp& tmp,
                                   double timeStepWidth,
                                   seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Volume integral
  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volume volKrnl = deviceVolumeKernelPrototype;
  kernel::gpu_localFlux localFluxKrnl = deviceLocalFluxKernelPrototype;

  const auto maxTmpMem = yateto::getMaxTmpMemRequired(volKrnl, localFluxKrnl);

  real* tmpMem = nullptr;
  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    unsigned maxNumElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    volKrnl.numElements = maxNumElements;

    // volume kernel always contains more elements than any local one
    tmpMem = (real*)(device.api->getStackMemory(maxTmpMem * maxNumElements));

    volKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    volKrnl.I =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    unsigned starOffset = 0;
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      volKrnl.star(i) =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      volKrnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }
    volKrnl.linearAllocator.initialize(tmpMem);
    volKrnl.streamPtr = runtime.stream();
    volKrnl.execute();
  }

  // Local Flux Integral
  for (unsigned face = 0; face < 4; ++face) {
    key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];
      localFluxKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      localFluxKrnl.I =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      localFluxKrnl.AplusT =
          const_cast<const real**>(entry.get(inner_keys::Wp::Id::AplusT)->getDeviceDataPtr());
      localFluxKrnl.linearAllocator.initialize(tmpMem);
      localFluxKrnl.streamPtr = runtime.stream();
      localFluxKrnl.execute(face);
    }

    ConditionalKey fsgKey(
        *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
    if (dataTable.find(fsgKey) != dataTable.end()) {
      auto nodalAvgDisplacements =
          dataTable[fsgKey].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();
      auto rhos = materialTable[fsgKey].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
      local_flux::aux::FreeSurfaceGravity freeSurfaceGravityBc;
      freeSurfaceGravityBc.g = gravitationalAcceleration;
      freeSurfaceGravityBc.rhos = rhos;
      freeSurfaceGravityBc.displacementDataPtrs = nodalAvgDisplacements;
      dirichletBoundary.evaluateOnDevice(face,
                                         fsgKey,
                                         deviceProjectRotatedKrnlPrototype,
                                         deviceNodalLfKrnlPrototype,
                                         freeSurfaceGravityBc,
                                         dataTable,
                                         device,
                                         runtime);
    }

    ConditionalKey dirichletKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
    if (dataTable.find(dirichletKey) != dataTable.end()) {
      auto easiBoundaryMapPtrs =
          dataTable[dirichletKey].get(inner_keys::Wp::Id::EasiBoundaryMap)->getDeviceDataPtr();
      auto easiBoundaryConstantPtrs =
          dataTable[dirichletKey].get(inner_keys::Wp::Id::EasiBoundaryConstant)->getDeviceDataPtr();

      local_flux::aux::EasiBoundary easiBoundaryBc;
      easiBoundaryBc.easiBoundaryMapPtrs = easiBoundaryMapPtrs;
      easiBoundaryBc.easiBoundaryConstantPtrs = easiBoundaryConstantPtrs;

      dirichletBoundary.evaluateOnDevice(face,
                                         dirichletKey,
                                         deviceProjectRotatedKrnlPrototype,
                                         deviceNodalLfKrnlPrototype,
                                         easiBoundaryBc,
                                         dataTable,
                                         device,
                                         runtime);
    }
  }
  if (tmpMem != nullptr) {
    device.api->popStackMemory();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Local::evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                           ConditionalIndicesTable& indicesTable,
                                           kernels::LocalData::Loader& loader,
                                           seissol::initializer::Layer& layer,
                                           seissol::initializer::LTS& lts,
                                           double time,
                                           double timeStepWidth,
                                           seissol::parallel::runtime::StreamRuntime& runtime) {

#ifdef ACL_DEVICE
  for (unsigned face = 0; face < 4; ++face) {
    ConditionalKey analyticalKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    if (indicesTable.find(analyticalKey) != indicesTable.end()) {
      const auto& cellIds =
          indicesTable[analyticalKey].get(inner_keys::Indices::Id::Cells)->getHostData();
      const size_t numElements = cellIds.size();
      auto* analytical = reinterpret_cast<real(*)[tensor::INodal::size()]>(
          layer.getScratchpadMemory(lts.analyticScratch));

      runtime.enqueueOmpFor(numElements, [=, &cellIds](std::size_t index) {
        auto cellId = cellIds.at(index);
        auto data = loader.entry(cellId);

        alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal::size()];

        assert(initConds != nullptr);
        assert(initConds->size() == 1);
        ApplyAnalyticalSolution applyAnalyticalSolution(this->getInitCond(0), data);

        dirichletBoundary.evaluateTimeDependent(nullptr,
                                                face,
                                                data.boundaryMapping()[face],
                                                m_projectKrnlPrototype,
                                                applyAnalyticalSolution,
                                                dofsFaceBoundaryNodal,
                                                time,
                                                timeStepWidth);

        std::memcpy(analytical[index], dofsFaceBoundaryNodal, sizeof(dofsFaceBoundaryNodal));
      });

      auto nodalLfKrnl = deviceNodalLfKrnlPrototype;
      nodalLfKrnl.INodal = const_cast<const real**>(
          dataTable[analyticalKey].get(inner_keys::Wp::Id::Analytical)->getDeviceDataPtr());
      nodalLfKrnl.AminusT = const_cast<const real**>(
          dataTable[analyticalKey].get(inner_keys::Wp::Id::AminusT)->getDeviceDataPtr());
      nodalLfKrnl.Q = dataTable[analyticalKey].get(inner_keys::Wp::Id::Dofs)->getDeviceDataPtr();
      nodalLfKrnl.streamPtr = runtime.stream();
      nodalLfKrnl.numElements = numElements;
      nodalLfKrnl.execute(face);
    }
  }
#else
  logError() << "No GPU implementation provided";
  ;
#endif // ACL_DEVICE
}

void Local::flopsIntegral(const FaceType faceTypes[4],
                          unsigned int& nonZeroFlops,
                          unsigned int& hardwareFlops) {
  nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for (unsigned int face = 0; face < 4; ++face) {
    // Local flux is executed for all faces that are not dynamic rupture.
    // For those cells, the flux is taken into account during the neighbor kernel.
    if (faceTypes[face] != FaceType::DynamicRupture) {
      nonZeroFlops += seissol::kernel::localFlux::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFlux::hardwareFlops(face);
    }

    // Take boundary condition flops into account.
    // Note that this only includes the flops of the kernels but not of the
    // boundary condition implementation.
    // The (probably incorrect) assumption is that they are negligible.
    switch (faceTypes[face]) {
    case FaceType::FreeSurfaceGravity:
      nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
                      seissol::kernel::projectToNodalBoundary::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
                       seissol::kernel::projectToNodalBoundary::hardwareFlops(face);
      break;
    case FaceType::Dirichlet:
      nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
                      seissol::kernel::projectToNodalBoundaryRotated::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
                       seissol::kernel::projectToNodalBoundary::hardwareFlops(face);
      break;
    case FaceType::Analytical:
      nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
                      ConvergenceOrder * seissol::kernel::updateINodal::NonZeroFlops;
      hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
                       ConvergenceOrder * seissol::kernel::updateINodal::HardwareFlops;
      break;
    default:
      break;
    }
  }
}

unsigned Local::bytesIntegral() {
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size();

  return reals * sizeof(real);
}

} // namespace seissol::kernels
