/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Local kernel of SeisSol.
 **/

#include "Kernels/Local.h"

#include <yateto.h>

#include <array>
#include <cassert>
#include <stdint.h>
#include "GravitationalFreeSurfaceBC.h"
#include "SeisSol.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop

#include <Kernels/common.hpp>
#include "../../../Kernels/Time.h"
GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

void seissol::kernels::LocalBase::checkGlobalData(GlobalData const* global, size_t alignment) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert(((uintptr_t)global->stiffnessMatrices(stiffness)) % alignment == 0);
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert(((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % alignment == 0);
    assert(((uintptr_t)global->changeOfBasisMatrices(flux)) % alignment == 0);
  }
#endif
}

void seissol::kernels::Local::setHostGlobalData(GlobalData const* global) {
  checkGlobalData(global, ALIGNMENT);
  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
  // for initial strain BC
  m_localInitFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localInitFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
  m_nodalLfKrnlPrototype.project2nFaceTo3m = global->project2nFaceTo3m;
  m_projectKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
  m_projectRotatedKrnlPrototype.V3mTo2nFace = global->V3mTo2nFace;
}

void seissol::kernels::Local::setGlobalData(const CompoundGlobalData& global) {
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
    initCondition->evaluate(time, nodesVec, localData.material, boundaryDofs);
  }

  private:
  seissol::physics::InitialField* initCondition{};
  LocalDataType& localData;
};

void seissol::kernels::Local::computeIntegral(
    real i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
    LocalData& data,
    LocalTmp& tmp,
    // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
    const CellMaterialData* materialData,
    CellBoundaryMapping const (*cellBoundaryMapping)[4],
    double time,
    double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(data.dofs) % ALIGNMENT == 0);

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.dofs;
  lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs + tensor::Q::size();

  for (int face = 0; face < 4; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
      lfKrnl.AplusT = data.localIntegration.nApNm1[face];
      if (data.cellInformation.faceTypes[face] != FaceType::regular &&
          data.cellInformation.faceTypes[face] != FaceType::periodic) {
        lfKrnl.execute(face);

        if (data.cellInformation.faceTypes[face] == FaceType::freeSurface ||
            data.cellInformation.faceTypes[face] == FaceType::outflow) {
          // additional term on free-surface BC to accomodate initial strain
          alignas(ALIGNMENT) real QInitialModal[tensor::Q::size()] = {0.0};
          alignas(ALIGNMENT) real QInitialNodal[tensor::QNodal::size()] = {0.0};
          real* exxNodal = (QInitialNodal + 0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          real* eyyNodal = (QInitialNodal + 1 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          real* ezzNodal = (QInitialNodal + 2 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          real* exyNodal = (QInitialNodal + 3 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          real* eyzNodal = (QInitialNodal + 4 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          real* ezxNodal = (QInitialNodal + 5 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
          for (unsigned int q = 0; q < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q) {
            // TODO(NONLINEAR) What are these numbers?
            exxNodal[q] = m_damagedElasticParameters->epsInitxx;
            eyyNodal[q] = m_damagedElasticParameters->epsInityy;
            ezzNodal[q] = m_damagedElasticParameters->epsInitzz;
            exyNodal[q] = m_damagedElasticParameters->epsInitxy;
            eyzNodal[q] = m_damagedElasticParameters->epsInityz;
            ezxNodal[q] = m_damagedElasticParameters->epsInitzx;
          }
          kernel::damageAssignFToDQ d_convertInitialToModal;
          d_convertInitialToModal.dQModal = QInitialModal;
          d_convertInitialToModal.vInv = init::vInv::Values;
          d_convertInitialToModal.FNodal = QInitialNodal;
          d_convertInitialToModal.execute();

          // Integrate it and add to Q
          kernel::localInitFlux lfIKrnl = m_localInitFluxKernelPrototype;
          lfIKrnl.Q = data.dofs;
          lfIKrnl.dQModal = QInitialModal;
          lfIKrnl.T = data.localIntegration.T[face];
          lfIKrnl.Tinv = data.localIntegration.Tinv[face];
          lfIKrnl.star(0) = data.localIntegration.ATtildeBC;
          lfIKrnl.fluxScale = data.localIntegration.fluxScales[face] * timeStepWidth;
          lfIKrnl.execute(face);
        }
      }
    }

    alignas(ALIGNMENT) real dofsFaceBoundaryNodal[tensor::INodal::size()];
    auto nodalLfKrnl = m_nodalLfKrnlPrototype;
    nodalLfKrnl.Q = data.dofs;
    nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
    nodalLfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
    nodalLfKrnl._prefetch.Q = data.dofs + tensor::Q::size();
    nodalLfKrnl.AminusT = data.neighboringIntegration.nAmNm1[face];

    // Include some boundary conditions here.
    switch (data.cellInformation.faceTypes[face]) {
    case FaceType::freeSurfaceGravity: {
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

      dirichletBoundary.evaluate(i_timeIntegratedDegreesOfFreedom,
                                 face,
                                 (*cellBoundaryMapping)[face],
                                 m_projectRotatedKrnlPrototype,
                                 applyFreeSurfaceBc,
                                 dofsFaceBoundaryNodal);

      nodalLfKrnl.execute(face);
      break;
    }
    case FaceType::dirichlet: {
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
      dirichletBoundary.evaluate(i_timeIntegratedDegreesOfFreedom,
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
    case FaceType::analytical: {
      assert(cellBoundaryMapping != nullptr);
      assert(initConds != nullptr);
      assert(initConds->size() == 1);
      ApplyAnalyticalSolution applyAnalyticalSolution(this->getInitCond(0), data);

      dirichletBoundary.evaluateTimeDependent(i_timeIntegratedDegreesOfFreedom,
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

void seissol::kernels::Local::computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                                     ConditionalMaterialTable& materialTable,
                                                     ConditionalIndicesTable& indicesTable,
                                                     kernels::LocalData::Loader& loader,
                                                     LocalTmp& tmp,
                                                     double timeStepWidth) {
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
    volKrnl.streamPtr = device.api->getDefaultStream();
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
      localFluxKrnl.streamPtr = device.api->getDefaultStream();
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
                                         device);
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
                                         device);
    }
  }
  if (tmpMem != nullptr) {
    device.api->popStackMemory();
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Local::evaluateBatchedTimeDependentBc(
    ConditionalPointersToRealsTable& dataTable,
    ConditionalIndicesTable& indicesTable,
    kernels::LocalData::Loader& loader,
    double time,
    double timeStepWidth) {

#ifdef ACL_DEVICE
  for (unsigned face = 0; face < 4; ++face) {
    ConditionalKey analyticalKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    if (indicesTable.find(analyticalKey) != indicesTable.end()) {
      auto idofsPtrs = dataTable[analyticalKey].get(inner_keys::Wp::Id::Idofs)->getHostData();

      auto cellIds = indicesTable[analyticalKey].get(inner_keys::Indices::Id::Cells)->getHostData();
      const size_t numElements = cellIds.size();

      for (unsigned index{0}; index < numElements; ++index) {
        auto cellId = cellIds[index];
        auto data = loader.entry(cellId);

        alignas(ALIGNMENT) real dofsFaceBoundaryNodal[tensor::INodal::size()];

        assert(initConds != nullptr);
        assert(initConds->size() == 1);
        ApplyAnalyticalSolution applyAnalyticalSolution(this->getInitCond(0), data);

        dirichletBoundary.evaluateTimeDependent(idofsPtrs[index],
                                                face,
                                                data.boundaryMapping[face],
                                                m_projectKrnlPrototype,
                                                applyAnalyticalSolution,
                                                dofsFaceBoundaryNodal,
                                                time,
                                                timeStepWidth);

        auto nodalLfKrnl = this->m_nodalLfKrnlPrototype;
        nodalLfKrnl.Q = data.dofs;
        nodalLfKrnl.INodal = dofsFaceBoundaryNodal;
        nodalLfKrnl.AminusT = data.neighboringIntegration.nAmNm1[face];
        nodalLfKrnl.execute(face);
      }
    }
  }
#else
  assert(false && "no implementation provided");
#endif // ACL_DEVICE
}

void seissol::kernels::Local::flopsIntegral(FaceType const i_faceTypes[4],
                                            unsigned int& o_nonZeroFlops,
                                            unsigned int& o_hardwareFlops) {
  o_nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  o_hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for (unsigned int face = 0; face < 4; ++face) {
    // Local flux is executed for all faces that are not dynamic rupture.
    // For those cells, the flux is taken into account during the neighbor kernel.
    if (i_faceTypes[face] != FaceType::dynamicRupture) {
      o_nonZeroFlops += seissol::kernel::localFlux::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFlux::hardwareFlops(face);
    }

    // Take boundary condition flops into account.
    // Note that this only includes the flops of the kernels but not of the
    // boundary condition implementation.
    // The (probably incorrect) assumption is that they are negligible.
    switch (i_faceTypes[face]) {
    case FaceType::freeSurfaceGravity:
      o_nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
                        seissol::kernel::projectToNodalBoundary::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
                         seissol::kernel::projectToNodalBoundary::hardwareFlops(face);
      break;
    case FaceType::dirichlet:
      o_nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
                        seissol::kernel::projectToNodalBoundaryRotated::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
                         seissol::kernel::projectToNodalBoundary::hardwareFlops(face);
      break;
    case FaceType::analytical:
      o_nonZeroFlops += seissol::kernel::localFluxNodal::nonZeroFlops(face) +
                        CONVERGENCE_ORDER * seissol::kernel::updateINodal::NonZeroFlops;
      o_hardwareFlops += seissol::kernel::localFluxNodal::hardwareFlops(face) +
                         CONVERGENCE_ORDER * seissol::kernel::updateINodal::HardwareFlops;
      break;
    default:
      break;
    }
  }
}

unsigned seissol::kernels::Local::bytesIntegral() {
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size();

  return reals * sizeof(real);
}

void seissol::kernels::Local::computeNonLinearRusanovFlux(
    const CellMaterialData* materialData,
    unsigned int l_cell,
    unsigned int side,
    const double* timeWeights,
    const real* qIPlus,
    const real* qIMinus,
    real* rusanovFluxP,
    const LocalIntegrationData* localIntegration) {
  using namespace seissol::dr::misc::quantity_indices;
  const unsigned DAM = 9;
  const unsigned BRE = 10;

  const real lambda0P = materialData[l_cell].local.lambda0;
  const real mu0P = materialData[l_cell].local.mu0;
  const real rho0P = materialData[l_cell].local.rho;

  const real lambda0M = materialData[l_cell].neighbor[side].lambda0;
  const real mu0M = materialData[l_cell].neighbor[side].mu0;
  const real rho0M = materialData[l_cell].neighbor[side].rho;

  const real epsInitxx = m_damagedElasticParameters->epsInitxx;
  const real epsInityy = m_damagedElasticParameters->epsInityy;
  const real epsInitzz = m_damagedElasticParameters->epsInitzz;
  const real epsInitxy = m_damagedElasticParameters->epsInitxy;
  const real epsInityz = m_damagedElasticParameters->epsInityz;
  const real epsInitzx = m_damagedElasticParameters->epsInitzx;

  const real aB0 = m_damagedElasticParameters->aB0;
  const real aB1 = m_damagedElasticParameters->aB1;
  const real aB2 = m_damagedElasticParameters->aB2;
  const real aB3 = m_damagedElasticParameters->aB3;

  real lambdaMax = 1.0 * std::sqrt((lambda0P + 2 * mu0P) / rho0P);
  real sxxP, syyP, szzP, sxyP, syzP, szxP, sxxM, syyM, szzM, sxyM, syzM, szxM;
  seissol::kernels::Time m_timeKernel;
  m_timeKernel.setDamagedElasticParameters(m_damagedElasticParameters);

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    auto weight = timeWeights[o];

    for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {

      real EspIp, EspIIp, xip, EspIm, EspIIm, xim;

      std::tie(EspIp, EspIIp, xip) = m_timeKernel.calculateEsp(
          &qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XX * seissol::dr::misc::numPaddedPoints],
          &qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  YY * seissol::dr::misc::numPaddedPoints],
          &qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  ZZ * seissol::dr::misc::numPaddedPoints],
          &qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XY * seissol::dr::misc::numPaddedPoints],
          &qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  YZ * seissol::dr::misc::numPaddedPoints],
          &qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XZ * seissol::dr::misc::numPaddedPoints],
          i,
          m_damagedElasticParameters);
      real alphap =
          qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                 DAM * seissol::dr::misc::numPaddedPoints + i];
      std::tie(EspIm, EspIIm, xim) = m_timeKernel.calculateEsp(
          &qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XX * seissol::dr::misc::numPaddedPoints],
          &qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   YY * seissol::dr::misc::numPaddedPoints],
          &qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   ZZ * seissol::dr::misc::numPaddedPoints],
          &qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XY * seissol::dr::misc::numPaddedPoints],
          &qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   YZ * seissol::dr::misc::numPaddedPoints],
          &qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XZ * seissol::dr::misc::numPaddedPoints],
          i,
          m_damagedElasticParameters);
      real alpham =
          qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  DAM * seissol::dr::misc::numPaddedPoints + i];
      real lambp =
          (1 - qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      BRE * seissol::dr::misc::numPaddedPoints + i]) *
              (lambda0P - alphap * materialData[l_cell].local.gammaR *
                              (qIPlus[o * seissol::dr::misc::numQuantities *
                                          seissol::dr::misc::numPaddedPoints +
                                      XX * seissol::dr::misc::numPaddedPoints + i] +
                               epsInitxx) /
                              std::sqrt(EspIIp)) +
          qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                 BRE * seissol::dr::misc::numPaddedPoints + i] *
              (2.0 * aB2 + 3.0 * xip * aB3 +
               aB1 *
                   (qIPlus[o * seissol::dr::misc::numQuantities *
                               seissol::dr::misc::numPaddedPoints +
                           XX * seissol::dr::misc::numPaddedPoints + i] +
                    epsInitxx) /
                   std::sqrt(EspIIp));
      real mup =
          (1 - qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      BRE * seissol::dr::misc::numPaddedPoints + i]) *
              (mu0P - alphap * materialData[l_cell].local.xi0 * materialData[l_cell].local.gammaR -
               0.5 * alphap * materialData[l_cell].local.gammaR * xip) +
          qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                 BRE * seissol::dr::misc::numPaddedPoints + i] *
              (aB0 + 0.5 * xip * aB1 - 0.5 * xip * xip * xip * aB3);

      real lambm =
          (1 - qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       BRE * seissol::dr::misc::numPaddedPoints + i]) *
              (lambda0M - alpham * materialData[l_cell].neighbor[side].gammaR *
                              (qIMinus[o * seissol::dr::misc::numQuantities *
                                           seissol::dr::misc::numPaddedPoints +
                                       XX * seissol::dr::misc::numPaddedPoints + i] +
                               epsInitxx) /
                              std::sqrt(EspIIm)) +
          qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  BRE * seissol::dr::misc::numPaddedPoints + i] *
              (2.0 * aB2 + 3.0 * xim * aB3 +
               aB1 *
                   (qIMinus[o * seissol::dr::misc::numQuantities *
                                seissol::dr::misc::numPaddedPoints +
                            XX * seissol::dr::misc::numPaddedPoints + i] +
                    epsInitxx) /
                   std::sqrt(EspIIm));

      real mum =
          (1 - qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       BRE * seissol::dr::misc::numPaddedPoints + i]) *
              (mu0M -
               alpham * materialData[l_cell].neighbor[side].xi0 *
                   materialData[l_cell].neighbor[side].gammaR -
               0.5 * alpham * materialData[l_cell].neighbor[side].gammaR * xim) +
          qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  BRE * seissol::dr::misc::numPaddedPoints + i] *
              (aB0 + 0.5 * xim * aB1 - 0.5 * xim * xim * xim * aB3);

      lambdaMax =
          std::min(std::sqrt((lambp + 2 * mup) / rho0P), std::sqrt((lambm + 2 * mum) / rho0M));

      // damage stress
      real mu_eff = materialData[l_cell].local.mu0 -
                    alphap * materialData[l_cell].local.gammaR * materialData[l_cell].local.xi0 -
                    0.5 * alphap * materialData[l_cell].local.gammaR * xip;
      real sxx_sp =
          materialData[l_cell].local.lambda0 * EspIp -
          alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
          2 * mu_eff *
              (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      XX * seissol::dr::misc::numPaddedPoints + i] +
               epsInitxx);
      real syy_sp =
          materialData[l_cell].local.lambda0 * EspIp -
          alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
          2 * mu_eff *
              (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      YY * seissol::dr::misc::numPaddedPoints + i] +
               epsInityy);
      real szz_sp =
          materialData[l_cell].local.lambda0 * EspIp -
          alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
          2 * mu_eff *
              (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      ZZ * seissol::dr::misc::numPaddedPoints + i] +
               epsInitzz);

      real sxy_sp =
          2 * mu_eff *
          (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XY * seissol::dr::misc::numPaddedPoints + i] +
           epsInitxy);
      real syz_sp =
          2 * mu_eff *
          (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  YZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInityz);
      real szx_sp =
          2 * mu_eff *
          (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInitzx);

      // breakage stress
      real sxx_bp =
          (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
          (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
              (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      XX * seissol::dr::misc::numPaddedPoints + i] +
               epsInitxx);
      real syy_bp =
          (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
          (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
              (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      YY * seissol::dr::misc::numPaddedPoints + i] +
               epsInityy);
      real szz_bp =
          (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
          (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
              (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                      ZZ * seissol::dr::misc::numPaddedPoints + i] +
               epsInitzz);

      real sxy_bp =
          (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
          (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XY * seissol::dr::misc::numPaddedPoints + i] +
           epsInitxy);
      real syz_bp =
          (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
          (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  YZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInityz);
      real szx_bp =
          (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) *
          (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  XZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInitzx);

      // damage stress minus
      mu_eff = materialData[l_cell].neighbor[side].mu0 -
               alpham * materialData[l_cell].neighbor[side].gammaR *
                   materialData[l_cell].neighbor[side].xi0 -
               0.5 * alpham * materialData[l_cell].neighbor[side].gammaR * xim;
      real sxx_sm =
          materialData[l_cell].neighbor[side].lambda0 * EspIm -
          alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
          2 * mu_eff *
              (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       XX * seissol::dr::misc::numPaddedPoints + i] +
               epsInitxx);
      real syy_sm =
          materialData[l_cell].neighbor[side].lambda0 * EspIm -
          alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
          2 * mu_eff *
              (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       YY * seissol::dr::misc::numPaddedPoints + i] +
               epsInityy);
      real szz_sm =
          materialData[l_cell].neighbor[side].lambda0 * EspIm -
          alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
          2 * mu_eff *
              (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       ZZ * seissol::dr::misc::numPaddedPoints + i] +
               epsInitzz);

      real sxy_sm =
          2 * mu_eff *
          (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XY * seissol::dr::misc::numPaddedPoints + i] +
           epsInitxy);
      real syz_sm =
          2 * mu_eff *
          (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   YZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInityz);
      real szx_sm =
          2 * mu_eff *
          (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInitzx);

      // breakage stress
      real sxx_bm =
          (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
          (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
              (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       XX * seissol::dr::misc::numPaddedPoints + i] +
               epsInitxx);
      real syy_bm =
          (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
          (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
              (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       YY * seissol::dr::misc::numPaddedPoints + i] +
               epsInityy);
      real szz_bm =
          (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
          (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
              (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       ZZ * seissol::dr::misc::numPaddedPoints + i] +
               epsInitzz);

      real sxy_bm =
          (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
          (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XY * seissol::dr::misc::numPaddedPoints + i] +
           epsInitxy);
      real syz_bm =
          (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
          (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   YZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInityz);
      real szx_bm =
          (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) *
          (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                   XZ * seissol::dr::misc::numPaddedPoints + i] +
           epsInitzx);

      real breakp =
          qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                 BRE * seissol::dr::misc::numPaddedPoints + i];
      real breakm =
          qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                  BRE * seissol::dr::misc::numPaddedPoints + i];

      sxxP = (1 - breakp) * sxx_sp + breakp * sxx_bp;
      syyP = (1 - breakp) * syy_sp + breakp * syy_bp;
      szzP = (1 - breakp) * szz_sp + breakp * szz_bp;
      sxyP = (1 - breakp) * sxy_sp + breakp * sxy_bp;
      syzP = (1 - breakp) * syz_sp + breakp * syz_bp;
      szxP = (1 - breakp) * szx_sp + breakp * szx_bp;

      sxxM = (1 - breakm) * sxx_sm + breakm * sxx_bm;
      syyM = (1 - breakm) * syy_sm + breakm * syy_bm;
      szzM = (1 - breakm) * szz_sm + breakm * szz_bm;

      sxyM = (1 - breakm) * sxy_sm + breakm * sxy_bm;
      syzM = (1 - breakm) * syz_sm + breakm * syz_bm;
      szxM = (1 - breakm) * szx_sm + breakm * szx_bm;

      rusanovFluxP[XX * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 *
                (-qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         U * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 * (-qIMinus[o * seissol::dr::misc::numQuantities *
                                seissol::dr::misc::numPaddedPoints +
                            U * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][0] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       XX * seissol::dr::misc::numPaddedPoints + i] +
                epsInitxx) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        XX * seissol::dr::misc::numPaddedPoints + i] +
                epsInitxx));

      rusanovFluxP[YY * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 *
                (-qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         V * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 * (-qIMinus[o * seissol::dr::misc::numQuantities *
                                seissol::dr::misc::numPaddedPoints +
                            V * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][1] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       YY * seissol::dr::misc::numPaddedPoints + i] +
                epsInityy) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        YY * seissol::dr::misc::numPaddedPoints + i] +
                epsInityy));

      rusanovFluxP[ZZ * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 *
                (-qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         W * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 * (-qIMinus[o * seissol::dr::misc::numQuantities *
                                seissol::dr::misc::numPaddedPoints +
                            W * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][2] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       ZZ * seissol::dr::misc::numPaddedPoints + i] +
                epsInitzz) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        ZZ * seissol::dr::misc::numPaddedPoints + i] +
                epsInitzz));

      rusanovFluxP[XY * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 *
                (-0.5 *
                 qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        V * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 *
                (-0.5 *
                 qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         V * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][0] +
           (0.5 *
                (-0.5 *
                 qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        U * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 *
                (-0.5 *
                 qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         U * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][1] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       XY * seissol::dr::misc::numPaddedPoints + i] +
                epsInitxy) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        XY * seissol::dr::misc::numPaddedPoints + i] +
                epsInitxy));

      rusanovFluxP[YZ * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 *
                (-0.5 *
                 qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        W * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 *
                (-0.5 *
                 qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         W * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][1] +
           (0.5 *
                (-0.5 *
                 qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        V * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 *
                (-0.5 *
                 qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         V * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][2] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       YZ * seissol::dr::misc::numPaddedPoints + i] +
                epsInityz) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        YZ * seissol::dr::misc::numPaddedPoints + i] +
                epsInityz));

      rusanovFluxP[XZ * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 *
                (-0.5 *
                 qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        W * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 *
                (-0.5 *
                 qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         W * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][0] +
           (0.5 *
                (-0.5 *
                 qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        U * seissol::dr::misc::numPaddedPoints + i]) +
            0.5 *
                (-0.5 *
                 qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                         U * seissol::dr::misc::numPaddedPoints + i])) *
               localIntegration[l_cell].surfaceNormal[side][2] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       XZ * seissol::dr::misc::numPaddedPoints + i] +
                epsInitzx) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        XZ * seissol::dr::misc::numPaddedPoints + i] +
                epsInitzx));

      rusanovFluxP[U * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 * (-sxxP / rho0P) + 0.5 * (-sxxM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][0] +
           (0.5 * (-sxyP / rho0P) + 0.5 * (-sxyM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][1] +
           (0.5 * (-szxP / rho0P) + 0.5 * (-szxM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][2] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       U * seissol::dr::misc::numPaddedPoints + i]) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        U * seissol::dr::misc::numPaddedPoints + i]));

      rusanovFluxP[V * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 * (-sxyP / rho0P) + 0.5 * (-sxyM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][0] +
           (0.5 * (-syyP / rho0P) + 0.5 * (-syyM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][1] +
           (0.5 * (-syzP / rho0P) + 0.5 * (-syzM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][2] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       V * seissol::dr::misc::numPaddedPoints + i]) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        V * seissol::dr::misc::numPaddedPoints + i]));

      rusanovFluxP[W * seissol::dr::misc::numPaddedPoints + i] +=
          weight *
          ((0.5 * (-szxP / rho0P) + 0.5 * (-szxM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][0] +
           (0.5 * (-syzP / rho0P) + 0.5 * (-syzM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][1] +
           (0.5 * (-szzP / rho0P) + 0.5 * (-szzM / rho0M)) *
               localIntegration[l_cell].surfaceNormal[side][2] +
           0.5 * lambdaMax *
               (qIPlus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                       W * seissol::dr::misc::numPaddedPoints + i]) -
           0.5 * lambdaMax *
               (qIMinus[o * seissol::dr::misc::numQuantities * seissol::dr::misc::numPaddedPoints +
                        W * seissol::dr::misc::numPaddedPoints + i]));
    }
  }
}

void seissol::kernels::Local::computeNonLinearIntegralCorrection(
    const CellLocalInformation* cellInformation,
    unsigned int l_cell,
    real** derivatives,
    real* (*faceNeighbors)[4],
    const CellMaterialData* materialData,
    const LocalIntegrationData* localIntegration,
    const NeighborData& data,
    const CellDRMapping (*drMapping)[4],
    kernel::nonlinearSurfaceIntegral& m_nonlSurfIntPrototype,
    double timeStepSize,
    const kernel::nonlEvaluateAndRotateQAtInterpolationPoints& m_nonlinearInterpolation) {
  double timePoints[CONVERGENCE_ORDER];
  double timeWeights[CONVERGENCE_ORDER];
  seissol::kernels::Time m_timeKernel;
  m_timeKernel.setDamagedElasticParameters(m_damagedElasticParameters);

  seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);

  for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
    timePoints[point] = 0.5 * (timeStepSize * timePoints[point] + timeStepSize);
    timeWeights[point] = 0.5 * timeStepSize * timeWeights[point];
  }
  for (unsigned int side = 0; side < 4; side++) {
    if (cellInformation[l_cell].faceTypes[side] == FaceType::regular ||
        cellInformation[l_cell].faceTypes[side] == FaceType::periodic) {
      // Compute local integrals with derivatives and Rusanov flux
      /// S1: compute the space-time interpolated Q on both side of 4 faces
      /// S2: at the same time rotate the field to face-aligned coord.
      alignas(PAGESIZE_STACK)
          real QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {
              {0.0}};
      alignas(PAGESIZE_STACK)
          real QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {
              {0.0}};

      for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
        real degreesOfFreedomPlus[tensor::Q::size()];
        real degreesOfFreedomMinus[tensor::Q::size()];

        for (unsigned i_f = 0; i_f < tensor::Q::size(); i_f++) {
          degreesOfFreedomPlus[i_f] = static_cast<real>(0.0);
          degreesOfFreedomMinus[i_f] = static_cast<real>(0.0);
        }

        // !!! Make sure every time after entering this function, the last input should be
        // reinitialized to zero
        m_timeKernel.computeTaylorExpansion(
            timePoints[timeInterval], 0.0, derivatives[l_cell], degreesOfFreedomPlus);
        m_timeKernel.computeTaylorExpansion(
            timePoints[timeInterval], 0.0, faceNeighbors[l_cell][side], degreesOfFreedomMinus);

        // Prototype is necessary for openmp
        kernel::nonlEvaluateAndRotateQAtInterpolationPoints m_nonLinInter =
            m_nonlinearInterpolation;

        m_nonLinInter.QInterpolated = &QInterpolatedPlus[timeInterval][0];
        m_nonLinInter.Q = degreesOfFreedomPlus;
        m_nonLinInter.execute(side, 0);

        m_nonLinInter.QInterpolated = &QInterpolatedMinus[timeInterval][0];
        m_nonLinInter.Q = degreesOfFreedomMinus;
        m_nonLinInter.execute(cellInformation[l_cell].faceRelations[side][0],
                              cellInformation[l_cell].faceRelations[side][1] + 1);
      }

      // S3: Construct matrices to store Rusanov flux on surface quadrature nodes.
      // Reshape the interpolated results
      using QInterpolatedShapeT =
          const real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];

      auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedPlus));
      auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedMinus));

      // The arrays to store time integrated flux
      alignas(PAGESIZE_STACK) real rusanovFluxPlus[tensor::QInterpolated::size()] = {0.0};
      // alignas(PAGESIZE_STACK) real rusanovFluxMinus[tensor::QInterpolated::size()] = {0.0};

      for (unsigned i_f = 0; i_f < tensor::QInterpolated::size(); i_f++) {
        rusanovFluxPlus[i_f] = static_cast<real>(0.0);
      }

      using rusanovFluxShape = real(*)[seissol::dr::misc::numPaddedPoints];
      auto* rusanovFluxP = reinterpret_cast<rusanovFluxShape>(rusanovFluxPlus);

      // S4: Compute the Rusanov flux
      computeNonLinearRusanovFlux(materialData,
                                  l_cell,
                                  side,
                                  timeWeights,
                                  *qIPlus[0],
                                  *qIMinus[0],
                                  *rusanovFluxP,
                                  localIntegration);

      /// S5: Integrate in space using quadrature.
      kernel::nonlinearSurfaceIntegral m_surfIntegral = m_nonlSurfIntPrototype;
      m_surfIntegral.Q = data.dofs;
      m_surfIntegral.Flux = rusanovFluxPlus;
      m_surfIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
      m_surfIntegral.execute(side, 0);
    } else if (cellInformation[l_cell].faceTypes[side] == FaceType::dynamicRupture) {
      // No neighboring cell contribution, interior bc.
      assert(reinterpret_cast<uintptr_t>(drMapping[l_cell][side].godunov) % ALIGNMENT == 0);

      kernel::nonlinearSurfaceIntegral m_drIntegral = m_nonlSurfIntPrototype;
      m_drIntegral.Q = data.dofs;
      m_drIntegral.Flux = drMapping[l_cell][side].godunov;
      m_drIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
      m_drIntegral.execute(side, drMapping[l_cell][side].faceRelation);
    } // if (faceTypes)
  }   // for (side)
}