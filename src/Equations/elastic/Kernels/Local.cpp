/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop

#include <Kernels/common.hpp>
GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

void seissol::kernels::LocalBase::checkGlobalData(GlobalData const* global, size_t alignment) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert( ((uintptr_t)global->stiffnessMatrices(stiffness)) % alignment == 0 );
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % alignment == 0 );
    assert( ((uintptr_t)global->changeOfBasisMatrices(flux)) % alignment == 0 );
  }
#endif
}

void seissol::kernels::Local::setHostGlobalData(GlobalData const* global) {
  checkGlobalData(global, ALIGNMENT);
  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;

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
  deviceLocalFluxKernelPrototype.plusFluxMatrices = global.onDevice->plusFluxMatrices;
  deviceNodalLfKrnlPrototype.project2nFaceTo3m = global.onDevice->project2nFaceTo3m;
  deviceProjectRotatedKrnlPrototype.V3mTo2nFace = global.onDevice->V3mTo2nFace;
#endif
}

template<typename LocalDataType>
struct ApplyAnalyticalSolution {
  ApplyAnalyticalSolution(seissol::physics::InitialField* initCondition,
                          LocalDataType& data) : initCondition(initCondition),
                                                 localData(data) {}


  void operator()(const real* nodes,
                  double time,
                  seissol::init::INodal::view::type& boundaryDofs) {

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

void seissol::kernels::Local::computeIntegral(real i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
                                              LocalData& data,
                                              LocalTmp& tmp,
                                              // TODO(Lukas) Nullable cause miniseissol. Maybe fix?
                                              const CellMaterialData* materialData,
                                              CellBoundaryMapping const (*cellBoundaryMapping)[4],
                                              double time,
                                              double timeStepWidth) {
  assert(reinterpret_cast<uintptr_t>(i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(data.dofs) % ALIGNMENT == 0);

  kernel::volume volKrnl = m_volumeKernelPrototype;
  volKrnl.Q = data.dofs;
  volKrnl.I = i_timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(volKrnl, get_ptr_sourceMatrix(data.localIntegration.specific));

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.dofs;
  lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs + tensor::Q::size();
  
  volKrnl.execute();

  for (int face = 0; face < 4; ++face) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if (data.cellInformation.faceTypes[face] != FaceType::dynamicRupture) {
      lfKrnl.AplusT = data.localIntegration.nApNm1[face];
      lfKrnl.execute(face);
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
    case FaceType::freeSurfaceGravity:
      {
        assert(cellBoundaryMapping != nullptr);
        assert(materialData != nullptr);
        auto* displ = tmp.nodalAvgDisplacements[face].data();
        auto displacement = init::averageNormalDisplacement::view::create(displ);
        auto applyFreeSurfaceBc = [&displacement, &materialData](
            const real*, // nodes are unused
            init::INodal::view::type& boundaryDofs) {
          for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
            const double rho = materialData->local.rho;
            const double g = getGravitationalAcceleration(); // [m/s^2]
            const double pressureAtBnd = -1 * rho * g * displacement(i);

            boundaryDofs(i,0) = 2 * pressureAtBnd - boundaryDofs(i,0);
            boundaryDofs(i,1) = 2 * pressureAtBnd - boundaryDofs(i,1);
            boundaryDofs(i,2) = 2 * pressureAtBnd - boundaryDofs(i,2);
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
    case FaceType::dirichlet:
      {
      assert(cellBoundaryMapping != nullptr);
      auto* easiBoundaryMap = (*cellBoundaryMapping)[face].easiBoundaryMap;
      auto* easiBoundaryConstant = (*cellBoundaryMapping)[face].easiBoundaryConstant;
      assert(easiBoundaryConstant != nullptr);
      assert(easiBoundaryMap != nullptr);
      auto applyEasiBoundary = [easiBoundaryMap, easiBoundaryConstant](
          const real* nodes,
          init::INodal::view::type& boundaryDofs) {
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
      case FaceType::analytical:
      {
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

void seissol::kernels::Local::computeBatchedIntegral(
  ConditionalPointersToRealsTable& dataTable,
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
    auto &entry = dataTable[key];

    unsigned maxNumElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    volKrnl.numElements = maxNumElements;

    // volume kernel always contains more elements than any local one
    tmpMem = (real*)(device.api->getStackMemory(maxTmpMem * maxNumElements));

    volKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    volKrnl.I = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());

    unsigned starOffset = 0;
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      volKrnl.star(i) = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
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
      auto &entry = dataTable[key];
      localFluxKrnl.numElements = entry.get(inner_keys::Wp::Id::Dofs)->getSize();
      localFluxKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      localFluxKrnl.I = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      localFluxKrnl.AplusT = const_cast<const real **>(entry.get(inner_keys::Wp::Id::AplusT)->getDeviceDataPtr());
      localFluxKrnl.linearAllocator.initialize(tmpMem);
      localFluxKrnl.streamPtr = device.api->getDefaultStream();
      localFluxKrnl.execute(face);
    }

    ConditionalKey fsgKey(*KernelNames::BoundaryConditions,
                          *ComputationKind::FreeSurfaceGravity,
                          face);
    if(dataTable.find(fsgKey) != dataTable.end()) {
      auto nodalAvgDisplacements = dataTable[fsgKey].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();
      auto rhos = materialTable[fsgKey].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
      local_flux::aux::FreeSurfaceGravity freeSurfaceGravityBc;
      freeSurfaceGravityBc.g = getGravitationalAcceleration();
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

    ConditionalKey dirichletKey(*KernelNames::BoundaryConditions,
                                *ComputationKind::Dirichlet,
                                face);
    if(dataTable.find(dirichletKey) != dataTable.end()) {
      auto easiBoundaryMapPtrs = dataTable[dirichletKey].get(inner_keys::Wp::Id::EasiBoundaryMap)->getDeviceDataPtr();
      auto easiBoundaryConstantPtrs = dataTable[dirichletKey].get(inner_keys::Wp::Id::EasiBoundaryConstant)->getDeviceDataPtr();

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
    ConditionalKey analyticalKey(*KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    if(indicesTable.find(analyticalKey) != indicesTable.end()) {
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
                                            unsigned int &o_nonZeroFlops,
                                            unsigned int &o_hardwareFlops)
{
  o_nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  o_hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for( unsigned int face = 0; face < 4; ++face ) {
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

unsigned seissol::kernels::Local::bytesIntegral()
{
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size();
  
  return reals * sizeof(real);
}
