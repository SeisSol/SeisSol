/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
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
 * Boundary kernel of SeisSol.
 **/

#include "Kernels/Neighbor.h"

#ifndef NDEBUG
#pragma message "compiling boundary kernel with assertions"
#endif

#include <cassert>
#include <stdint.h>

void seissol::kernels::NeighborBase::checkGlobalData(GlobalData const* global, size_t alignment) {
#ifndef NDEBUG
  for( int l_neighbor = 0; l_neighbor < 4; ++l_neighbor ) {
    assert( ((uintptr_t)global->changeOfBasisMatrices(l_neighbor)) % alignment == 0 );
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(l_neighbor)) % alignment == 0 );
    assert( ((uintptr_t)global->neighbourChangeOfBasisMatricesTransposed(l_neighbor)) % alignment == 0 );
  }

  for( int h = 0; h < 3; ++h ) {
    assert( ((uintptr_t)global->neighbourFluxMatrices(h)) % alignment == 0 );
  }

  for (int i = 0; i < 4; ++i) {
    for(int h = 0; h < 3; ++h) {
      assert( ((uintptr_t)global->nodalFluxMatrices(i,h)) % alignment == 0 );
    }
  }
#endif
}

void seissol::kernels::Neighbor::setHostGlobalData(GlobalData const* global) {
  checkGlobalData(global, ALIGNMENT);
  m_nfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global->neighbourChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global->neighbourFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global->nodalFluxMatrices;
}

void seissol::kernels::Neighbor::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);

  deviceNfKrnlPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceNfKrnlPrototype.rT = global.onDevice->neighbourChangeOfBasisMatricesTransposed;
  deviceNfKrnlPrototype.fP = global.onDevice->neighbourFluxMatrices;
  deviceDrKrnlPrototype.V3mTo2nTWDivM = global.onDevice->nodalFluxMatrices;
#endif
}

void seissol::kernels::Neighbor::computeNeighborsIntegral(NeighborData& data,
                                                          CellDRMapping const (&cellDrMapping)[4],
                                                          real* i_timeIntegrated[4],
                                                          real* faceNeighbors_prefetch[4]) {
  assert(reinterpret_cast<uintptr_t>(data.dofs) % ALIGNMENT == 0);

  for (unsigned int l_face = 0; l_face < 4; l_face++) {
    switch (data.cellInformation.faceTypes[l_face]) {
    case FaceType::regular:
      // Fallthrough intended
    case FaceType::periodic:
      {
      // Standard neighboring flux
      // Compute the neighboring elements flux matrix id.
      assert(reinterpret_cast<uintptr_t>(i_timeIntegrated[l_face]) % ALIGNMENT == 0 );
      assert(data.cellInformation.faceRelations[l_face][0] < 4
             && data.cellInformation.faceRelations[l_face][1] < 3);
      kernel::neighboringFlux nfKrnl = m_nfKrnlPrototype;
      nfKrnl.Q = data.dofs;
      nfKrnl.I = i_timeIntegrated[l_face];
      nfKrnl.AminusT = data.neighboringIntegration.nAmNm1[l_face];
      nfKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
      nfKrnl.execute(data.cellInformation.faceRelations[l_face][1],
		     data.cellInformation.faceRelations[l_face][0],
		     l_face);
      break;
      }
    case FaceType::dynamicRupture:
      {
      // No neighboring cell contribution, interior bc.
      assert(reinterpret_cast<uintptr_t>(cellDrMapping[l_face].godunov) % ALIGNMENT == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[l_face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[l_face].godunov;
      drKrnl.Q = data.dofs;
      drKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
      drKrnl.execute(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
      break;
      }
    default:
      // No contribution for all other cases.
      // Note: some other bcs are handled in the local kernel.
      break;
    }
  }
}

void seissol::kernels::Neighbor::computeBatchedNeighborsIntegral(ConditionalBatchTableT &table) {
#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype;

  real* tmpMem = nullptr;
  device.api->resetCircularStreamCounter();
  auto resetDeviceCurrentState = [this](size_t counter) {
    for (size_t i = 0; i < counter; ++i) {
      this->device.api->popStackMemory();
    }
    this->device.api->fastStreamsSync();
    this->device.api->resetCircularStreamCounter();
  };


  for(size_t face = 0; face < 4; face++) {
    size_t streamCounter{0};

    // regular and periodic
    for (size_t faceRelation = 0; faceRelation < (*FaceRelations::Count); ++faceRelation) {

      ConditionalKey key(*KernelNames::NeighborFlux,
                         (FaceKinds::Regular || FaceKinds::Periodic),
                         face,
                         faceRelation);

      if(table.find(key) != table.end()) {
        BatchTable &entry = table[key];

        const auto NUM_ELEMENTS = (entry.content[*EntityId::Dofs])->getSize();
        neighFluxKrnl.numElements = NUM_ELEMENTS;

        neighFluxKrnl.Q = (entry.content[*EntityId::Dofs])->getPointers();
        neighFluxKrnl.I = const_cast<const real **>((entry.content[*EntityId::Idofs])->getPointers());
        neighFluxKrnl.AminusT = const_cast<const real **>((entry.content[*EntityId::AminusT])->getPointers());
        neighFluxKrnl.streamPtr = device.api->getNextCircularStream();

        tmpMem = (real*)(device.api->getStackMemory(neighFluxKrnl.TmpMaxMemRequiredInBytes * NUM_ELEMENTS));
        neighFluxKrnl.linearAllocator.initialize(tmpMem);

        (neighFluxKrnl.*neighFluxKrnl.ExecutePtrs[faceRelation])();
        ++streamCounter;
      }
    }

    // dynamic rupture
    for (unsigned faceRelation = 0; faceRelation < (*DrFaceRelations::Count); ++faceRelation) {

      ConditionalKey Key(*KernelNames::NeighborFlux,
                         *FaceKinds::DynamicRupture,
                         face,
                         faceRelation);

      if(table.find(Key) != table.end()) {
        BatchTable &entry = table[Key];

        const auto NUM_ELEMENTS = (entry.content[*EntityId::Dofs])->getSize();
        drKrnl.numElements = NUM_ELEMENTS;

        drKrnl.fluxSolver = const_cast<const real **>((entry.content[*EntityId::FluxSolver])->getPointers());
        drKrnl.QInterpolated = const_cast<real const**>((entry.content[*EntityId::Godunov])->getPointers());
        drKrnl.Q = (entry.content[*EntityId::Dofs])->getPointers();
        drKrnl.streamPtr = device.api->getNextCircularStream();

        tmpMem = (real*)(device.api->getStackMemory(drKrnl.TmpMaxMemRequiredInBytes * NUM_ELEMENTS));
        drKrnl.linearAllocator.initialize(tmpMem);

        (drKrnl.*drKrnl.ExecutePtrs[faceRelation])();
        ++streamCounter;
      }
    }
    resetDeviceCurrentState(streamCounter);
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Neighbor::flopsNeighborsIntegral(const FaceType i_faceTypes[4],
                                                        const int i_neighboringIndices[4][2],
                                                        CellDRMapping const (&cellDrMapping)[4],
                                                        unsigned int &o_nonZeroFlops,
                                                        unsigned int &o_hardwareFlops,
                                                        long long& o_drNonZeroFlops,
                                                        long long& o_drHardwareFlops) {
  // reset flops
  o_nonZeroFlops = 0;
  o_hardwareFlops = 0;
  o_drNonZeroFlops = 0;
  o_drHardwareFlops = 0;
  
  for (unsigned int face = 0; face < 4; face++) {
    // compute the neighboring elements flux matrix id.
    switch (i_faceTypes[face]) {
    case FaceType::regular:
      // Fallthrough intended
    case FaceType::periodic:
      // regular neighbor
      assert(i_neighboringIndices[face][0] < 4 && i_neighboringIndices[face][1] < 3);
      o_nonZeroFlops += kernel::neighboringFlux::nonZeroFlops(i_neighboringIndices[face][1], i_neighboringIndices[face][0], face);
      o_hardwareFlops += kernel::neighboringFlux::hardwareFlops(i_neighboringIndices[face][1], i_neighboringIndices[face][0], face);
      break;
    case FaceType::dynamicRupture:
      o_drNonZeroFlops += dynamicRupture::kernel::nodalFlux::nonZeroFlops(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      o_drHardwareFlops += dynamicRupture::kernel::nodalFlux::hardwareFlops(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      break;
    default:
      //Handled in local kernel
      break;
    }
  }
}


unsigned seissol::kernels::Neighbor::bytesNeighborsIntegral()
{
  unsigned reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size();
  
  return reals * sizeof(real);
}
