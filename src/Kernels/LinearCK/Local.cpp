// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Kernels/LinearCK/LocalBase.h"

#include <Alignment.h>
#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <GeneratedCode/init.h>
#include <GeneratedCode/kernel.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Parallel/Runtime/Stream.h>
#include <Physics/InitialField.h>
#include <Solver/MultipleSimulations.h>
#include <cstddef>
#include <cstdint>
#include <memory>
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

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)
namespace seissol::kernels::solver::linearck {

  template<typename Cfg>
void Local<Cfg>::setGlobalData(const CompoundGlobalData& global) {
  m_volumeKernelPrototype.kDivM = global.onHost->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global.onHost->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global.onHost->localChangeOfBasisMatricesTransposed;

  m_nodalLfKrnlPrototype.project2nFaceTo3m = global.onHost->project2nFaceTo3m;

  m_projectKrnlPrototype.V3mTo2nFace = global.onHost->v3mTo2nFace;
  m_projectRotatedKrnlPrototype.V3mTo2nFace = global.onHost->v3mTo2nFace;

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();

  deviceVolumeKernelPrototype.kDivM = global.onDevice->stiffnessMatrices;
#ifdef USE_PREMULTIPLY_FLUX
  deviceLocalFluxKernelPrototype.plusFluxMatrices = global.onDevice->plusFluxMatrices;
#else
  deviceLocalFluxKernelPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceLocalFluxKernelPrototype.fMrT = global.onDevice->localChangeOfBasisMatricesTransposed;
#endif
  deviceNodalLfKrnlPrototype.project2nFaceTo3m = global.onDevice->project2nFaceTo3m;
  deviceProjectRotatedKrnlPrototype.V3mTo2nFace = global.onDevice->v3mTo2nFace;
#endif
}

template<typename Cfg>
struct ApplyAnalyticalSolution {
  ApplyAnalyticalSolution(const std::vector<std::unique_ptr<physics::InitialField>>* initConditions,
                          LTS::Ref<Cfg>& data)
      : initConditions(initConditions), localData(data) {}

  void operator()(const real* nodes, double time, typename seissol::init::INodal<Cfg>::view::type& boundaryDofs) {
    assert(initConditions != nullptr);

    constexpr auto NodeCount = seissol::tensor::INodal<Cfg>::Shape[multisim::BasisFunctionDimension];
    alignas(Alignment) std::array<double, 3> nodesVec[NodeCount];

#pragma omp simd
    for (std::size_t i = 0; i < NodeCount; ++i) {
      nodesVec[i][0] = nodes[i * 3 + 0];
      nodesVec[i][1] = nodes[i * 3 + 1];
      nodesVec[i][2] = nodes[i * 3 + 2];
    }

    for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
      auto slicedBoundaryDofs = multisim::simtensor(boundaryDofs, s);
      initConditions->at(s % initConditions->size())
          ->evaluate(time, nodesVec, NodeCount, localData.template get<LTS::Material>(), slicedBoundaryDofs);
    }
  }

  private:
  const std::vector<std::unique_ptr<physics::InitialField>>* initConditions;
  LTS::Ref<Cfg>& localData;
};

template<typename Cfg>
void Local<Cfg>::computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I<Cfg>::size()],
                            LTS::Ref<Cfg>& data,
                            LocalTmp<Cfg>& tmp,
                            // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                            const CellMaterialData* materialData,
                            const CellBoundaryMapping (*cellBoundaryMapping)[4],
                            double time,
                            double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(timeIntegratedDegreesOfFreedom) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>()) % Alignment == 0);

  kernel::volume<Cfg> volKrnl = m_volumeKernelPrototype;
  volKrnl.Q = data.template get<LTS::Dofs>();
  volKrnl.I = timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
    volKrnl.star(i) = data.template get<LTS::LocalIntegration>().starMatrices[i];
  }

  // Optional source term
  set_ET(volKrnl, get_ptr_sourceMatrix(data.template get<LTS::LocalIntegration>().specific));

  kernel::localFlux<Cfg> lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.template get<LTS::Dofs>();
  lfKrnl.I = timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I<Cfg>::size();
  lfKrnl._prefetch.Q = data.template get<LTS::Dofs>() + tensor::Q<Cfg>::size();

  volKrnl.execute();

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.template get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
      lfKrnl.AplusT = data.template get<LTS::LocalIntegration>().nApNm1[face];
      lfKrnl.execute(face);
    }

    alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal<Cfg>::size()];
    auto nodalLfKrnl = m_nodalLfKrnlPrototype;
    nodalLfKrnl.Q = data.template get<LTS::Dofs>();
    nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
    nodalLfKrnl._prefetch.I = timeIntegratedDegreesOfFreedom + tensor::I<Cfg>::size();
    nodalLfKrnl._prefetch.Q = data.template get<LTS::Dofs>() + tensor::Q<Cfg>::size();
    nodalLfKrnl.AminusT = data.template get<LTS::NeighboringIntegration>().nAmNm1[face];

    // Include some boundary conditions here.
    switch (data.template get<LTS::CellInformation>().faceTypes[face]) {
    case FaceType::FreeSurfaceGravity: {
      assert(cellBoundaryMapping != nullptr);
      assert(materialData != nullptr);
      auto* displ = tmp.nodalAvgDisplacements[face].data();
      auto displacement = init::averageNormalDisplacement<Cfg>::view::create(displ);
      // lambdas can't catch gravitationalAcceleration directly, so have to make a copy here.
      const auto localG = this->gravitationalAcceleration;
      auto applyFreeSurfaceBc =
          [&](const real*, // nodes are unused
              typename init::INodal<Cfg>::view::type& boundaryDofs) {
            for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
              auto slicedBoundaryDofs = multisim::simtensor(boundaryDofs, s);
              auto slicedDisplacement = multisim::simtensor(displacement, s);

              for (unsigned int i = 0;
                   i < nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisFunctionDimension];
                   ++i) {
                const double rho = materialData->local->getDensity();
                assert(localG > 0);
                const double pressureAtBnd = -1 * rho * localG * slicedDisplacement(i);

                slicedBoundaryDofs(i, 0) = 2 * pressureAtBnd - slicedBoundaryDofs(i, 0);
                slicedBoundaryDofs(i, 1) = 2 * pressureAtBnd - slicedBoundaryDofs(i, 1);
                slicedBoundaryDofs(i, 2) = 2 * pressureAtBnd - slicedBoundaryDofs(i, 2);
              }
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
                                   const real* nodes, typename init::INodal<Cfg>::view::type& boundaryDofs) {
        seissol::kernel::createEasiBoundaryGhostCells<Cfg> easiBoundaryKernel;
        easiBoundaryKernel.easiBoundaryMap = easiBoundaryMap;
        easiBoundaryKernel.easiBoundaryConstant = easiBoundaryConstant;
        easiBoundaryKernel.easiIdentMap = init::easiIdentMap<Cfg>::Values;
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
      assert(this->initConds != nullptr);
      ApplyAnalyticalSolution<Cfg> applyAnalyticalSolution(this->initConds, data);

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

template<typename Cfg>
void Local<Cfg>::computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                   ConditionalMaterialTable& materialTable,
                                   ConditionalIndicesTable& indicesTable,
                                   double timeStepWidth,
                                   seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Volume integral
  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volume<Cfg> volKrnl = deviceVolumeKernelPrototype;
  kernel::gpu_localFlux<Cfg> localFluxKrnl = deviceLocalFluxKernelPrototype;

  const auto maxTmpMem = yateto::getMaxTmpMemRequired(volKrnl, localFluxKrnl);

  // volume kernel always contains more elements than any local one
  const auto maxNumElements = dataTable.find(key) != dataTable.end()
                                  ? (dataTable[key].get(inner_keys::Wp::Id::Dofs))->getSize()
                                  : 0;
  auto tmpMem = runtime.memoryHandle<real>((maxTmpMem * maxNumElements) / sizeof(real));
  if (dataTable.find(key) != dataTable.end()) {
    auto& entry = dataTable[key];

    volKrnl.numElements = maxNumElements;

    volKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    volKrnl.I =
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      volKrnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      volKrnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }
    volKrnl.linearAllocator.initialize(tmpMem.get());
    volKrnl.streamPtr = runtime.stream();
    volKrnl.execute();
  }

  // Local Flux Integral
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];
      localFluxKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      localFluxKrnl.I =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      localFluxKrnl.AplusT = const_cast<const real**>(
          entry.get(inner_keys::Wp::Id::LocalIntegrationData)->getDeviceDataPtr());
      localFluxKrnl.extraOffset_AplusT = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, nApNm1, face);
      localFluxKrnl.linearAllocator.initialize(tmpMem.get());
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
#else
  logError() << "No GPU implementation provided";
#endif
}

template<typename Cfg>
void Local<Cfg>::evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                           ConditionalIndicesTable& indicesTable,
                                           LTS::Layer& layer,
                                           double time,
                                           double timeStepWidth,
                                           seissol::parallel::runtime::StreamRuntime& runtime) {

#ifdef ACL_DEVICE
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    ConditionalKey analyticalKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    if (indicesTable.find(analyticalKey) != indicesTable.end()) {
      const auto& cellIds =
          indicesTable[analyticalKey].get(inner_keys::Indices::Id::Cells)->getHostData();
      const size_t numElements = cellIds.size();
      auto* analytical =
          reinterpret_cast<real(*)[tensor::INodal<Cfg>::size()]>(layer.var<LTS::AnalyticScratch>());

      runtime.enqueueLoop(numElements, [=, &cellIds, &layer](std::size_t index) {
        auto cellId = cellIds.at(index);
        auto data = layer.cellRef<Cfg>(cellId);

        alignas(Alignment) real dofsFaceBoundaryNodal[tensor::INodal<Cfg>::size()];

        assert(this->initConds != nullptr);
        ApplyAnalyticalSolution<Cfg> applyAnalyticalSolution(this->initConds, data);

        dirichletBoundary.evaluateTimeDependent(nullptr,
                                                face,
                                                data.template get<LTS::BoundaryMapping>()[face],
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
      nodalLfKrnl.AminusT =
          const_cast<const real**>(dataTable[analyticalKey]
                                       .get(inner_keys::Wp::Id::NeighborIntegrationData)
                                       ->getDeviceDataPtr());
      nodalLfKrnl.extraOffset_AminusT =
          SEISSOL_ARRAY_OFFSET(NeighboringIntegrationData, nAmNm1, face);
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

template<typename Cfg>
void Local<Cfg>::flopsIntegral(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                          std::uint64_t& nonZeroFlops,
                          std::uint64_t& hardwareFlops) {
  nonZeroFlops = seissol::kernel::volume<Cfg>::NonZeroFlops;
  hardwareFlops = seissol::kernel::volume<Cfg>::HardwareFlops;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // Local flux is executed for all faces that are not dynamic rupture.
    // For those cells, the flux is taken into account during the neighbor kernel.
    if (faceTypes[face] != FaceType::DynamicRupture) {
      nonZeroFlops += seissol::kernel::localFlux<Cfg>::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFlux<Cfg>::hardwareFlops(face);
    }

    // Take boundary condition flops into account.
    // Note that this only includes the flops of the kernels but not of the
    // boundary condition implementation.
    // The (probably incorrect) assumption is that they are negligible.
    switch (faceTypes[face]) {
    case FaceType::FreeSurfaceGravity:
      nonZeroFlops += seissol::kernel::localFluxNodal<Cfg>::nonZeroFlops(face) +
                      seissol::kernel::projectToNodalBoundary<Cfg>::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxNodal<Cfg>::hardwareFlops(face) +
                       seissol::kernel::projectToNodalBoundary<Cfg>::hardwareFlops(face);
      break;
    case FaceType::Dirichlet:
      nonZeroFlops += seissol::kernel::localFluxNodal<Cfg>::nonZeroFlops(face) +
                      seissol::kernel::projectToNodalBoundaryRotated<Cfg>::nonZeroFlops(face);
      hardwareFlops += seissol::kernel::localFluxNodal<Cfg>::hardwareFlops(face) +
                       seissol::kernel::projectToNodalBoundary<Cfg>::hardwareFlops(face);
      break;
    case FaceType::Analytical:
      nonZeroFlops += seissol::kernel::localFluxNodal<Cfg>::nonZeroFlops(face) +
                      Cfg::ConvergenceOrder * seissol::kernel::updateINodal<Cfg>::NonZeroFlops;
      hardwareFlops += seissol::kernel::localFluxNodal<Cfg>::hardwareFlops(face) +
                       Cfg::ConvergenceOrder * seissol::kernel::updateINodal<Cfg>::HardwareFlops;
      break;
    default:
      break;
    }
  }
}

template<typename Cfg>
std::uint64_t Local<Cfg>::bytesIntegral() {
  std::uint64_t reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star<Cfg>>();
  // flux solvers
  reals += static_cast<std::uint64_t>(4 * tensor::AplusT<Cfg>::size());

  // DOFs write
  reals += tensor::Q<Cfg>::size();

  return reals * sizeof(real);
}

#define _H_(cfg) template class Local<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels::solver::linearck
