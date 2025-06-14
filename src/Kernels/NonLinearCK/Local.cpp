// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "Kernels/NonLinearCK/LocalBase.h"

// Required for material update
#include <Equations/Setup.h>

// #include "LocalBase.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>

// Required for material update
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include <Parallel/Runtime/Stream.h>
#include <Physics/InitialField.h>
#include <Solver/MultipleSimulations.h>
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

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)
namespace seissol::kernels::solver::nonlinearck {

void setStarMatrix(
    const real* matAT, const real* matBT, const real* matCT, const real grad[3], real* starMatrix) {
  for (unsigned idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
    starMatrix[idx] = grad[0] * matAT[idx];
  }

  for (unsigned idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
    starMatrix[idx] += grad[1] * matBT[idx];
  }

  for (unsigned idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
    starMatrix[idx] += grad[2] * matCT[idx];
  }
}

void Local::setGlobalData(const CompoundGlobalData& global) {
  m_volumeKernelPrototype.kDivM = global.onHost->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global.onHost->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global.onHost->localChangeOfBasisMatricesTransposed;

  m_nodalLfKrnlPrototype.project2nFaceTo3m = global.onHost->project2nFaceTo3m;

  m_projectKrnlPrototype.V3mTo2nFace = global.onHost->V3mTo2nFace;
  m_projectRotatedKrnlPrototype.V3mTo2nFace = global.onHost->V3mTo2nFace;

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
    for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
      auto slicedBoundaryDofs = multisim::simtensor(boundaryDofs, s);

      for (unsigned int i = 0; i < seissol::tensor::INodal::Shape[0]; ++i) {
        auto curNode = std::array<double, 3>{};
        curNode[0] = nodes[offset++];
        curNode[1] = nodes[offset++];
        curNode[2] = nodes[offset++];
        nodesVec.push_back(curNode);
      }

      assert(initCondition != nullptr);
      initCondition->evaluate(time, nodesVec, localData.material(), slicedBoundaryDofs);
    }
  }

  private:
  seissol::physics::InitialField* initCondition{};
  LocalDataType& localData;
};

void Local::updateMaterials(LocalData& data) {
  // Summary of the requried input:
  // LocalData: material.local, material.neighbor, cellInformation, dofs,
  // LocalData.specific: localVertices, localElements, meshId,
  // Update local modulus
  // logWarning() << data.cellInformation().faceTypes[0] << "\n";
  // kernels::LocalData::Loader loader;
  // loader.load(*m_lts, i_layerData);

  // CellMaterialData* materialData = i_layerData.var(m_lts->material);

  real ATData[tensor::star::size(0)];
  real ATtildeData[tensor::star::size(0)];
  real BTData[tensor::star::size(1)];
  real CTData[tensor::star::size(2)];
  auto AT = init::star::view<0>::create(ATData);
  // AT with elastic parameters in local coordinate system, used for flux kernel
  auto ATtilde = init::star::view<0>::create(ATtildeData);
  auto BT = init::star::view<0>::create(BTData);
  auto CT = init::star::view<0>::create(CTData);

  real TData[seissol::tensor::T::size()];
  real TinvData[seissol::tensor::Tinv::size()];
  auto T = init::T::view::create(TData);
  auto Tinv = init::Tinv::view::create(TinvData);

  real QgodLocalData[tensor::QgodLocal::size()];
  real QgodNeighborData[tensor::QgodNeighbor::size()];
  auto QgodLocal = init::QgodLocal::view::create(QgodLocalData);
  auto QgodNeighbor = init::QgodNeighbor::view::create(QgodNeighborData);

  real rusanovPlusNull[tensor::QcorrLocal::size()]{};
  real rusanovMinusNull[tensor::QcorrNeighbor::size()]{};

  // kernel::cellAve m_cellAverageKernel;
  // // real Q_aveData[NUMBER_OF_QUANTITIES];
  // real Q_aveData[tensor::QAve::size()];
  // auto Q_ave = init::QAve::view::create(Q_aveData);>

  // data = loader.entry(l_cell);

  // // TODO: compute cell-average of solutions and compute new material properties
  // // real Q_ave[NUMBER_OF_QUANTITIES];
  // m_cellAverageKernel.phiAve = init::phiAve::Values;
  // m_cellAverageKernel.Q = data.dofs;
  // m_cellAverageKernel.QAve = Q_aveData;
  // m_cellAverageKernel.execute();

  // Change to updated material properties
  data.material().local.mu = data.material().local.muE * 0.9;
  data.material().local.lambda = data.material().local.lambdaE * 0.9;
  for (unsigned i_s = 0; i_s < 4; i_s++) {
    data.material().neighbor[i_s].mu = data.material().neighbor[i_s].muE * 0.9;
    data.material().neighbor[i_s].lambda = data.material().neighbor[i_s].lambdaE * 0.9;
  }

  /// global coordinates of the vertices
  real x[4];
  real y[4];
  real z[4];
  real gradXi[3];
  real gradEta[3];
  real gradZeta[3];

  // Iterate over all 4 vertices of the tetrahedron
  for (unsigned vertex = 0; vertex < 4; ++vertex) {
    const VrtxCoords& coords = data.localIntegration().specific.localVertices[vertex].coords;
    x[vertex] = coords[0];
    y[vertex] = coords[1];
    z[vertex] = coords[2];
  }

  seissol::transformations::tetrahedronGlobalToReferenceJacobian(
      x, y, z, gradXi, gradEta, gradZeta);
  seissol::model::getTransposedCoefficientMatrix(data.material().local, 0, AT);
  seissol::model::getTransposedCoefficientMatrix(data.material().local, 1, BT);
  seissol::model::getTransposedCoefficientMatrix(data.material().local, 2, CT);
  setStarMatrix(ATData, BTData, CTData, gradXi, data.localIntegration().starMatrices[0]);
  setStarMatrix(ATData, BTData, CTData, gradEta, data.localIntegration().starMatrices[1]);
  setStarMatrix(ATData, BTData, CTData, gradZeta, data.localIntegration().starMatrices[2]);

  real volume = data.localIntegration().specific.localVolume;

  for (unsigned side = 0; side < 4; ++side) {
    VrtxCoords normal;
    VrtxCoords tangent1;
    VrtxCoords tangent2;
    for (unsigned coord = 0; coord < 3; ++coord) {
      normal[coord] = data.localIntegration().specific.localNormal[side][coord];
      tangent1[coord] = data.localIntegration().specific.localTangent1[side][coord];
      tangent2[coord] = data.localIntegration().specific.localTangent2[side][coord];
    }
    double surface = data.localIntegration().specific.localSurfaces[side];

    // real NLocalData[6*6];
    // seissol::model::getBondMatrix(normal, tangent1, tangent2, NLocalData);
    // if (materialData[l_cell].local.getMaterialType() ==
    // seissol::model::MaterialType::anisotropic) {
    //   seissol::model::getTransposedGodunovState(
    //   seissol::model::getRotatedMaterialCoefficients(NLocalData,
    //   *dynamic_cast<seissol::model::AnisotropicMaterial*>(&materialData[l_cell].local)),
    //                                               seissol::model::getRotatedMaterialCoefficients(NLocalData,
    //                                               *dynamic_cast<seissol::model::AnisotropicMaterial*>(&materialData[l_cell].neighbor[side])),
    //                                               cellInformation[l_cell].faceTypes[side],
    //                                               QgodLocal,
    //                                               QgodNeighbor );
    //   seissol::model::getTransposedCoefficientMatrix(
    //   seissol::model::getRotatedMaterialCoefficients(NLocalData,
    //   *dynamic_cast<seissol::model::AnisotropicMaterial*>(&materialData[l_cell].local)), 0,
    //   ATtilde );
    // } else {
    seissol::model::getTransposedGodunovState(data.material().local,
                                              data.material().neighbor[side],
                                              data.cellInformation().faceTypes[side],
                                              QgodLocal,
                                              QgodNeighbor);
    seissol::model::getTransposedCoefficientMatrix(data.material().local, 0, ATtilde);
    // }

    // Calculate transposed T instead
    seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, T, Tinv);

    // Scale with |S_side|/|J| and multiply with -1 as the flux matrices
    // must be subtracted.
    real fluxScale = -2.0 * surface / (6.0 * volume);

    kernel::computeFluxSolverLocal localKrnl;
    localKrnl.fluxScale = fluxScale;
    localKrnl.AplusT = data.localIntegration().nApNm1[side];
    localKrnl.QgodLocal = QgodLocalData;
    localKrnl.QcorrLocal = rusanovPlusNull;
    localKrnl.T = TData;
    localKrnl.Tinv = TinvData;
    localKrnl.star(0) = ATtildeData;
    localKrnl.execute();

    kernel::computeFluxSolverNeighbor neighKrnl;
    neighKrnl.fluxScale = fluxScale;
    neighKrnl.AminusT = data.neighboringIntegration().nAmNm1[side];
    neighKrnl.QgodNeighbor = QgodNeighborData;
    neighKrnl.QcorrNeighbor = rusanovMinusNull;
    neighKrnl.T = TData;
    neighKrnl.Tinv = TinvData;
    neighKrnl.star(0) = ATtildeData;
    if (data.cellInformation().faceTypes[side] == FaceType::Dirichlet ||
        data.cellInformation().faceTypes[side] == FaceType::FreeSurfaceGravity) {
      // Already rotated!
      neighKrnl.Tinv = init::identityT::Values;
    }
    neighKrnl.execute();
  }
}

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
            for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
              auto slicedBoundaryDofs = multisim::simtensor(boundaryDofs, s);
              auto slicedDisplacement = multisim::simtensor(displacement, s);

              for (unsigned int i = 0;
                   i < nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension];
                   ++i) {
                const double rho = materialData->local.rho;
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
#ifdef USE_DAMAGE
  assert(reinterpret_cast<uintptr_t>(timeIntegratedDegreesOfFreedom) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(data.dofs()) % Alignment == 0);
  // Do volumetric flux integration (volKrnl) + local face flux integration (lfKrnl)
  // One additional: nonlinear source term
  // Previous TimeCluster.cpp
#endif
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

    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      volKrnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      volKrnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }
    volKrnl.linearAllocator.initialize(tmpMem.get());
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

} // namespace seissol::kernels::solver::nonlinearck
