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
 * Neighbor kernel of SeisSol.
 **/

#include "Kernels/Neighbor.h"

#include <cassert>
#include <stdint.h>
#include <cstddef>
#include <cstring>

#include <generated_code/init.h>

void seissol::kernels::Neighbor::setHostGlobalData(GlobalData const* global) {
#ifndef NDEBUG
  for( int l_neighbor = 0; l_neighbor < 4; ++l_neighbor ) {
    assert( ((uintptr_t)global->changeOfBasisMatrices(l_neighbor)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(l_neighbor)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->neighborChangeOfBasisMatricesTransposed(l_neighbor)) % ALIGNMENT == 0 );
  }

  for( int h = 0; h < 3; ++h ) {
    assert( ((uintptr_t)global->neighborFluxMatrices(h)) % ALIGNMENT == 0 );
  }

  for (int i = 0; i < 4; ++i) {
    for(int h = 0; h < 3; ++h) {
      assert( ((uintptr_t)global->nodalFluxMatrices(i,h)) % ALIGNMENT == 0 );
    }
  }
#endif
  m_nfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global->neighborChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global->neighborFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global->nodalFluxMatrices;
  m_nKrnlPrototype.selectEla = init::selectEla::Values;
  m_nKrnlPrototype.selectAne = init::selectAne::Values;
}

void seissol::kernels::Neighbor::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);
}

void seissol::kernels::Neighbor::computeNeighborsIntegral(  NeighborData&                     data,
                                                            CellDRMapping const             (&cellDrMapping)[4],
                                                            real*                             i_timeIntegrated[4],
                                                            real*                             faceNeighbors_prefetch[4] )
{
#ifndef NDEBUG
  for( int l_neighbor = 0; l_neighbor < 4; ++l_neighbor ) {
    // alignment of the time integrated dofs
    if( data.cellInformation().faceTypes[l_neighbor] != FaceType::outflow && data.cellInformation().faceTypes[l_neighbor] != FaceType::dynamicRupture ) { // no alignment for outflow and DR boundaries required
      assert( ((uintptr_t)i_timeIntegrated[l_neighbor]) % ALIGNMENT == 0 );
    }
  }
#endif

  // alignment of the degrees of freedom
  assert( ((uintptr_t)data.dofs()) % ALIGNMENT == 0 );

  real Qext[tensor::Qext::size()] __attribute__((aligned(PAGESIZE_STACK))) = {};

  kernel::neighborFluxExt nfKrnl = m_nfKrnlPrototype;
  nfKrnl.Qext = Qext;

  // iterate over faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( data.cellInformation().faceTypes[l_face] != FaceType::outflow
        && data.cellInformation().faceTypes[l_face] != FaceType::dynamicRupture ) {
      // compute the neighboring elements flux matrix id.
      if( data.cellInformation().faceTypes[l_face] != FaceType::freeSurface ) {
        assert(data.cellInformation().faceRelations[l_face][0] < 4 && data.cellInformation().faceRelations[l_face][1] < 3);

        nfKrnl.I = i_timeIntegrated[l_face];
        nfKrnl.AminusT = data.neighboringIntegration().nAmNm1[l_face];
        nfKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
        nfKrnl.execute(data.cellInformation().faceRelations[l_face][1], data.cellInformation().faceRelations[l_face][0], l_face);
      }
    } else if (data.cellInformation().faceTypes[l_face] == FaceType::dynamicRupture) {
      assert(((uintptr_t)cellDrMapping[l_face].godunov) % ALIGNMENT == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[l_face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[l_face].godunov;
      drKrnl.Qext = Qext;
      drKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
      drKrnl.execute(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
    }
  }

  kernel::neighbor nKrnl = m_nKrnlPrototype;
  nKrnl.Qext = Qext;
  nKrnl.Q = data.dofs();
  nKrnl.Qane = data.dofsAne();
  nKrnl.w = data.neighboringIntegration().specific.w;

  nKrnl.execute();
}

void seissol::kernels::Neighbor::flopsNeighborsIntegral(const FaceType i_faceTypes[4],
                                                        const int i_neighboringIndices[4][2],
                                                        CellDRMapping const (&cellDrMapping)[4],
                                                        unsigned int &o_nonZeroFlops,
                                                        unsigned int &o_hardwareFlops,
                                                        long long& o_drNonZeroFlops,
                                                        long long& o_drHardwareFlops) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;
  o_drNonZeroFlops = 0; o_drHardwareFlops = 0;

  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if(i_faceTypes[l_face] != FaceType::outflow
       && i_faceTypes[l_face] != FaceType::dynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if( i_faceTypes[l_face] != FaceType::freeSurface ) {
        assert(i_neighboringIndices[l_face][0] < 4 && i_neighboringIndices[l_face][1] < 3);

        o_nonZeroFlops  += seissol::kernel::neighborFluxExt::nonZeroFlops(i_neighboringIndices[l_face][1], i_neighboringIndices[l_face][0], l_face);
        o_hardwareFlops += seissol::kernel::neighborFluxExt::hardwareFlops(i_neighboringIndices[l_face][1], i_neighboringIndices[l_face][0], l_face);
      } else { // fall back to local matrices in case of free surface boundary conditions
        o_nonZeroFlops  += seissol::kernel::localFluxExt::nonZeroFlops(l_face);
        o_hardwareFlops += seissol::kernel::localFluxExt::hardwareFlops(l_face);
      }
    } else if (i_faceTypes[l_face] == FaceType::dynamicRupture) {
      o_drNonZeroFlops += dynamicRupture::kernel::nodalFlux::nonZeroFlops(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
      o_drHardwareFlops += dynamicRupture::kernel::nodalFlux::hardwareFlops(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
    }
  }

  o_nonZeroFlops += kernel::neighbor::NonZeroFlops;
  o_hardwareFlops += kernel::neighbor::HardwareFlops;
}


unsigned seissol::kernels::Neighbor::bytesNeighborsIntegral()
{
  unsigned reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size() + 2 * tensor::Qane::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size() + tensor::w::size();

  return reals * sizeof(real);
}

void seissol::kernels::Neighbor::computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable &table) {
#ifdef ACL_DEVICE
  kernel::gpu_neighboringFluxExt neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype;

  real* tmpMem = nullptr;
  device.api->resetCircularStreamCounter();
  auto resetDeviceCurrentState = [this](size_t counter) {
    for (size_t i = 0; i < counter; ++i) {
      this->device.api->popStackMemory();
    }
    this->device.api->joinCircularStreamsToDefault();
    this->device.api->resetCircularStreamCounter();
  };

  for(size_t face = 0; face < 4; face++) {
    this->device.api->forkCircularStreamsFromDefault();
    size_t streamCounter{0};

    // regular and periodic
    for (size_t faceRelation = 0; faceRelation < (*FaceRelations::Count); ++faceRelation) {

      ConditionalKey key(*KernelNames::NeighborFlux,
                         (FaceKinds::Regular || FaceKinds::Periodic),
                         face,
                         faceRelation);

      if(table.find(key) != table.end()) {
        auto &entry = table[key];

        const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
        neighFluxKrnl.numElements = numElements;

        neighFluxKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
        neighFluxKrnl.I = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
        neighFluxKrnl.AminusT = const_cast<const real **>((entry.get(inner_keys::Wp::Id::AminusT))->getDeviceDataPtr());

        tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(neighFluxKrnl.TmpMaxMemRequiredInBytes * numElements));
        neighFluxKrnl.linearAllocator.initialize(tmpMem);

        neighFluxKrnl.streamPtr = device.api->getNextCircularStream();
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
        auto &entry = table[Key];

        const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
        drKrnl.numElements = numElements;

        drKrnl.fluxSolver = const_cast<const real **>((entry.get(inner_keys::Wp::Id::FluxSolver))->getDeviceDataPtr());
        drKrnl.QInterpolated = const_cast<real const**>((entry.get(inner_keys::Wp::Id::Godunov))->getDeviceDataPtr());
        drKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

        tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(drKrnl.TmpMaxMemRequiredInBytes * numElements));
        drKrnl.linearAllocator.initialize(tmpMem);

        drKrnl.streamPtr = device.api->getNextCircularStream();
        (drKrnl.*drKrnl.ExecutePtrs[faceRelation])();
        ++streamCounter;
      }
    }
    resetDeviceCurrentState(streamCounter);

    ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    if (dataTable.find(key) != dataTable.end()) {
      kernel::gpu_neighbor nKrnl;
      nKrnl.numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
      nKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();
      nKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      nKrnl.Qane = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
      nKrnl.w = (entry.get(inner_keys::Wp::Id::Omega))->getDeviceDataPtr();

      nKrnl.execute();
    }
  }
#else
  assert(false && "no implementation provided");
#endif
}
