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

#include "generated_code/init.h"

namespace seissol::kernels {

void Neighbor::setHostGlobalData(GlobalData const* global) {
#ifndef NDEBUG
  for( int neighbor = 0; neighbor < 4; ++neighbor ) {
    assert( (reinterpret_cast<uintptr_t>(global->changeOfBasisMatrices(neighbor))) % Alignment == 0 );
    assert( (reinterpret_cast<uintptr_t>(global->localChangeOfBasisMatricesTransposed(neighbor))) % Alignment == 0 );
    assert( (reinterpret_cast<uintptr_t>(global->neighborChangeOfBasisMatricesTransposed(neighbor))) % Alignment == 0 );
  }

  for( int h = 0; h < 3; ++h ) {
    assert( (reinterpret_cast<uintptr_t>(global->neighborFluxMatrices(h))) % Alignment == 0 );
  }

  for (int i = 0; i < 4; ++i) {
    for(int h = 0; h < 3; ++h) {
      assert( (reinterpret_cast<uintptr_t>(global->nodalFluxMatrices(i,h))) % Alignment == 0 );
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

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);
#ifdef ACL_DEVICE
#ifdef USE_PREMULTIPLY_FLUX
  deviceNfKrnlPrototype.minusFluxMatrices = global.onDevice->minusFluxMatrices;
#else
  deviceNfKrnlPrototype.rDivM = global.onDevice->changeOfBasisMatrices;
  deviceNfKrnlPrototype.rT = global.onDevice->neighborChangeOfBasisMatricesTransposed;
  deviceNfKrnlPrototype.fP = global.onDevice->neighborFluxMatrices;
#endif
  deviceDrKrnlPrototype.V3mTo2nTWDivM = global.onDevice->nodalFluxMatrices;
  deviceNKrnlPrototype.selectEla = global.onDevice->selectEla;
  deviceNKrnlPrototype.selectAne = global.onDevice->selectAne;
#endif
}

void Neighbor::computeNeighborsIntegral(  NeighborData&                     data,
                                                            CellDRMapping const             (&cellDrMapping)[4],
                                                            real*                             timeIntegrated[4],
                                                            real*                             faceNeighborsPrefetch[4] )
{
#ifndef NDEBUG
  for( int neighbor = 0; neighbor < 4; ++neighbor ) {
    // alignment of the time integrated dofs
    if( data.cellInformation().faceTypes[neighbor] != FaceType::Outflow && data.cellInformation().faceTypes[neighbor] != FaceType::DynamicRupture ) { // no alignment for outflow and DR boundaries required
      assert( (reinterpret_cast<uintptr_t>(timeIntegrated[neighbor])) % Alignment == 0 );
    }
  }
#endif

  // alignment of the degrees of freedom
  assert( (reinterpret_cast<uintptr_t>(data.dofs())) % Alignment == 0 );

  alignas(PagesizeStack) real qext[tensor::Qext::size()] = {};

  kernel::neighborFluxExt nfKrnl = m_nfKrnlPrototype;
  nfKrnl.Qext = qext;

  // iterate over faces
  for( unsigned int face = 0; face < 4; face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( data.cellInformation().faceTypes[face] != FaceType::Outflow
        && data.cellInformation().faceTypes[face] != FaceType::DynamicRupture ) {
      // compute the neighboring elements flux matrix id.
      if( data.cellInformation().faceTypes[face] != FaceType::FreeSurface ) {
        assert(data.cellInformation().faceRelations[face][0] < 4 && data.cellInformation().faceRelations[face][1] < 3);

        nfKrnl.I = timeIntegrated[face];
        nfKrnl.AminusT = data.neighboringIntegration().nAmNm1[face];
        nfKrnl._prefetch.I = faceNeighborsPrefetch[face];
        nfKrnl.execute(data.cellInformation().faceRelations[face][1], data.cellInformation().faceRelations[face][0], face);
      }
    } else if (data.cellInformation().faceTypes[face] == FaceType::DynamicRupture) {
      assert((reinterpret_cast<uintptr_t>(cellDrMapping[face].godunov)) % Alignment == 0);

      dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[face].fluxSolver;
      drKrnl.QInterpolated = cellDrMapping[face].godunov;
      drKrnl.Qext = qext;
      drKrnl._prefetch.I = faceNeighborsPrefetch[face];
      drKrnl.execute(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  kernel::neighbor nKrnl = m_nKrnlPrototype;
  nKrnl.Qext = qext;
  nKrnl.Q = data.dofs();
  nKrnl.Qane = data.dofsAne();
  nKrnl.w = data.neighboringIntegration().specific.w;

  nKrnl.execute();
}

void Neighbor::flopsNeighborsIntegral(const FaceType faceTypes[4],
                                                        const int neighboringIndices[4][2],
                                                        CellDRMapping const (&cellDrMapping)[4],
                                                        unsigned int &nonZeroFlops,
                                                        unsigned int &hardwareFlops,
                                                        long long& drNonZeroFlops,
                                                        long long& drHardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;
  drNonZeroFlops = 0; drHardwareFlops = 0;

  for( unsigned int face = 0; face < 4; face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if(faceTypes[face] != FaceType::Outflow
       && faceTypes[face] != FaceType::DynamicRupture) {
      // compute the neighboring elements flux matrix id.
      if( faceTypes[face] != FaceType::FreeSurface ) {
        assert(neighboringIndices[face][0] < 4 && neighboringIndices[face][1] < 3);

        nonZeroFlops  += seissol::kernel::neighborFluxExt::nonZeroFlops(neighboringIndices[face][1], neighboringIndices[face][0], face);
        hardwareFlops += seissol::kernel::neighborFluxExt::hardwareFlops(neighboringIndices[face][1], neighboringIndices[face][0], face);
      } else { // fall back to local matrices in case of free surface boundary conditions
        nonZeroFlops  += seissol::kernel::localFluxExt::nonZeroFlops(face);
        hardwareFlops += seissol::kernel::localFluxExt::hardwareFlops(face);
      }
    } else if (faceTypes[face] == FaceType::DynamicRupture) {
      drNonZeroFlops += dynamicRupture::kernel::nodalFlux::nonZeroFlops(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
      drHardwareFlops += dynamicRupture::kernel::nodalFlux::hardwareFlops(cellDrMapping[face].side, cellDrMapping[face].faceRelation);
    }
  }

  nonZeroFlops += kernel::neighbor::NonZeroFlops;
  hardwareFlops += kernel::neighbor::HardwareFlops;
}


unsigned Neighbor::bytesNeighborsIntegral()
{
  unsigned reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * tensor::I::size() + 2 * tensor::Q::size() + 2 * tensor::Qane::size();
  // flux solvers load
  reals += 4 * tensor::AminusT::size() + tensor::w::size();

  return reals * sizeof(real);
}

void Neighbor::computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable &table, seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt neighFluxKrnl = deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux drKrnl = deviceDrKrnlPrototype;

  {
    ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    if (table.find(key) != table.end()) {
      auto &entry = table[key];
      device.algorithms.setToValue((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr(),
                                        static_cast<real>(0.0),
                                        tensor::Qext::Size,
                                        (entry.get(inner_keys::Wp::Id::DofsExt))->getSize(),
                                        runtime.stream());
    }
  }

  real* tmpMem = nullptr;
  auto resetDeviceCurrentState = [this](size_t counter) {
    for (size_t i = 0; i < counter; ++i) {
      this->device.api->popStackMemory();
    }
  };

  for(size_t face = 0; face < 4; face++) {
    std::size_t streamCounter = 0;
    runtime.envMany((*FaceRelations::Count)+(*DrFaceRelations::Count), [&](void* stream, size_t i) {
      // regular and periodic
      if (i < (*FaceRelations::Count)) {
        // regular and periodic
        unsigned faceRelation = i;

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

          tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(seissol::kernel::gpu_neighborFluxExt::TmpMaxMemRequiredInBytes * numElements));
          neighFluxKrnl.linearAllocator.initialize(tmpMem);

          neighFluxKrnl.streamPtr = stream;
          (neighFluxKrnl.*seissol::kernel::gpu_neighborFluxExt::ExecutePtrs[faceRelation])();
          ++streamCounter;
        }
      }
      else {
        // Dynamic Rupture
        unsigned faceRelation = i - (*FaceRelations::Count);

        ConditionalKey key(*KernelNames::NeighborFlux,
                          *FaceKinds::DynamicRupture,
                          face,
                          faceRelation);

        if(table.find(key) != table.end()) {
          auto &entry = table[key];

          const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
          drKrnl.numElements = numElements;

          drKrnl.fluxSolver = const_cast<const real **>((entry.get(inner_keys::Wp::Id::FluxSolver))->getDeviceDataPtr());
          drKrnl.QInterpolated = const_cast<real const**>((entry.get(inner_keys::Wp::Id::Godunov))->getDeviceDataPtr());
          drKrnl.Qext = (entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr();

          tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(seissol::dynamicRupture::kernel::gpu_nodalFlux::TmpMaxMemRequiredInBytes * numElements));
          drKrnl.linearAllocator.initialize(tmpMem);

          drKrnl.streamPtr = stream;
          (drKrnl.*seissol::dynamicRupture::kernel::gpu_nodalFlux::ExecutePtrs[faceRelation])();
          ++streamCounter;
        }
      }
    });
    resetDeviceCurrentState(streamCounter);
  }

  ConditionalKey key(KernelNames::Time || KernelNames::Volume);
  if (table.find(key) != table.end()) {
    auto &entry = table[key];
    kernel::gpu_neighbor nKrnl = deviceNKrnlPrototype;
    nKrnl.numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    nKrnl.Qext = const_cast<const real **>((entry.get(inner_keys::Wp::Id::DofsExt))->getDeviceDataPtr());
    nKrnl.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
    nKrnl.Qane = (entry.get(inner_keys::Wp::Id::DofsAne))->getDeviceDataPtr();
    nKrnl.w = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Omega))->getDeviceDataPtr());
    nKrnl.streamPtr = runtime.stream();

    nKrnl.execute();
  }
#else
  assert(false && "no implementation provided");
#endif
}
} // namespace seissol::kernels

