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

#ifndef NDEBUG
#pragma message "compiling local kernel with assertions"
#endif

#include <yateto.h>


#include <array>
#include <cassert>
#include <stdint.h>

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
  deviceLocalFluxKernelPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceLocalFluxKernelPrototype.fMrT = global.onDevice->localChangeOfBasisMatricesTransposed;
  deviceNodalLfKrnlPrototype.project2nFaceTo3m = global.onDevice->project2nFaceTo3m;
#endif
}

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
            const double g = 9.81; // [m/s^2]
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
      auto applyAnalyticalSolution = [materialData, this](const real* nodes,
                                                    double time,
                                                    init::INodal::view::type& boundaryDofs) {
          auto nodesVec = std::vector<std::array<double, 3>>{};
          int offset = 0;
          for (unsigned int i = 0; i < tensor::INodal::Shape[0]; ++i) {
            auto curNode = std::array<double, 3>{};
            curNode[0] = nodes[offset++];
            curNode[1] = nodes[offset++];
            curNode[2] = nodes[offset++];
            nodesVec.push_back(curNode);
          }
          assert(initConds != nullptr);
          // TODO(Lukas) Support multiple init. conds?
          assert(initConds->size() == 1);
          (*initConds)[0]->evaluate(time, nodesVec, *materialData, boundaryDofs);
      };

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

void seissol::kernels::Local::computeBatchedIntegral(ConditionalBatchTableT &table, LocalTmp& tmp) {
#ifdef ACL_DEVICE
  // Volume integral
  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  kernel::gpu_volume volKrnl = deviceVolumeKernelPrototype;
  kernel::gpu_localFlux localFluxKrnl = deviceLocalFluxKernelPrototype;

  constexpr size_t MAX_TMP_MEM = (volKrnl.TmpMaxMemRequiredInBytes > localFluxKrnl.TmpMaxMemRequiredInBytes) \
                                   ? volKrnl.TmpMaxMemRequiredInBytes : localFluxKrnl.TmpMaxMemRequiredInBytes;

  real* tmpMem = nullptr;
  if (table.find(key) != table.end()) {
    BatchTable &entry = table[key];

    unsigned MaxNumElements = (entry.content[*EntityId::Dofs])->getSize();
    volKrnl.numElements = MaxNumElements;

    // volume kernel always contains more elements than any local one
    tmpMem = (real*)(device.api->getStackMemory(MAX_TMP_MEM * MaxNumElements));

    volKrnl.Q = (entry.content[*EntityId::Dofs])->getPointers();
    volKrnl.I = const_cast<const real **>((entry.content[*EntityId::Idofs])->getPointers());

    unsigned starOffset = 0;
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      volKrnl.star(i) = const_cast<const real **>((entry.content[*EntityId::Star])->getPointers());
      volKrnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }
    volKrnl.linearAllocator.initialize(tmpMem);
    volKrnl.execute();
  }

  // Local Flux Integral
  for (unsigned face = 0; face < 4; ++face) {
    key = ConditionalKey(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);

    if (table.find(key) != table.end()) {
      BatchTable &entry = table[key];
      localFluxKrnl.numElements = entry.content[*EntityId::Dofs]->getSize();
      localFluxKrnl.Q = (entry.content[*EntityId::Dofs])->getPointers();
      localFluxKrnl.I = const_cast<const real **>((entry.content[*EntityId::Idofs])->getPointers());
      localFluxKrnl.AplusT = const_cast<const real **>(entry.content[*EntityId::AplusT]->getPointers());
      localFluxKrnl.linearAllocator.initialize(tmpMem);
      localFluxKrnl.execute(face);
    }
  }
  if (tmpMem != nullptr) {
    device.api->popStackMemory();
  }
#else
  assert(false && "no implementation provided");
#endif
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
