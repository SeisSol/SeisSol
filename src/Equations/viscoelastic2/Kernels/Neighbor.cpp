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

#ifndef NDEBUG
#pragma message "compiling boundary kernel with assertions"
#endif

#include <cassert>
#include <stdint.h>
#include <cstddef>
#include <cstring>

#include <generated_code/init.h>

void seissol::kernels::Neighbor::setGlobalData(GlobalData const* global) {
#ifndef NDEBUG
  for( int l_neighbor = 0; l_neighbor < 4; ++l_neighbor ) {
    assert( ((uintptr_t)global->changeOfBasisMatrices(l_neighbor)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(l_neighbor)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->neighbourChangeOfBasisMatricesTransposed(l_neighbor)) % ALIGNMENT == 0 );
  }
  
  for( int h = 0; h < 3; ++h ) {
    assert( ((uintptr_t)global->neighbourFluxMatrices(h)) % ALIGNMENT == 0 );
  }
  
  for (int i = 0; i < 4; ++i) {
    for(int h = 0; h < 3; ++h) {
      assert( ((uintptr_t)global->nodalFluxMatrices(i,h)) % ALIGNMENT == 0 );
    }
  }
#endif
  m_lfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
  m_lfKrnlPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.rDivM = global->changeOfBasisMatrices;
  m_nfKrnlPrototype.rT = global->neighbourChangeOfBasisMatricesTransposed;
  m_nfKrnlPrototype.fP = global->neighbourFluxMatrices;
  m_drKrnlPrototype.V3mTo2nTWDivM = global->nodalFluxMatrices;
  m_nKrnlPrototype.selectEla = init::selectEla::Values;
  m_nKrnlPrototype.selectAne = init::selectAne::Values;
}

void seissol::kernels::Neighbor::computeNeighborsIntegral(  NeighborData&                     data,
                                                            CellDRMapping const             (&cellDrMapping)[4],
                                                            real*                             i_timeIntegrated[4],
                                                            real*                             faceNeighbors_prefetch[4] )
{
#ifndef NDEBUG
  for( int l_neighbor = 0; l_neighbor < 4; ++l_neighbor ) {
    // alignment of the time integrated dofs
    if( data.cellInformation.faceTypes[l_neighbor] != outflow && data.cellInformation.faceTypes[l_neighbor] != dynamicRupture ) { // no alignment for outflow and DR boundaries required
      assert( ((uintptr_t)i_timeIntegrated[l_neighbor]) % ALIGNMENT == 0 );
    }
  }
#endif

  // alignment of the degrees of freedom
  assert( ((uintptr_t)data.dofs) % ALIGNMENT == 0 );

  real Qext[tensor::Qext::size()] __attribute__((aligned(PAGESIZE_STACK))) = {};

  kernel::neighbourFluxExt nfKrnl = m_nfKrnlPrototype;
  nfKrnl.Qext = Qext;

  // iterate over faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( data.cellInformation.faceTypes[l_face] != outflow && data.cellInformation.faceTypes[l_face] != dynamicRupture ) {
      // compute the neighboring elements flux matrix id.
      if( data.cellInformation.faceTypes[l_face] != freeSurface ) {
        assert(data.cellInformation.faceRelations[l_face][0] < 4 && data.cellInformation.faceRelations[l_face][1] < 3);
        
        nfKrnl.I = i_timeIntegrated[l_face];
        nfKrnl.AminusT = data.neighboringIntegration.nAmNm1[l_face];
        nfKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
        nfKrnl.execute(data.cellInformation.faceRelations[l_face][1], data.cellInformation.faceRelations[l_face][0], l_face);
      } else { // fall back to local matrices in case of free surface boundary conditions
        kernel::localFluxExt lfKrnl = m_lfKrnlPrototype;
        lfKrnl.Qext = data.dofs;
        lfKrnl.I = i_timeIntegrated[l_face];
        lfKrnl.AplusT = data.neighboringIntegration.nAmNm1[l_face];
        lfKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
        lfKrnl.execute(l_face);
      }
    } else if (data.cellInformation.faceTypes[l_face] == dynamicRupture) {
      assert(((uintptr_t)cellDrMapping[l_face].godunov) % ALIGNMENT == 0);

      kernel::nodalFlux drKrnl = m_drKrnlPrototype;
      drKrnl.fluxSolver = cellDrMapping[l_face].fluxSolver;
      drKrnl.godunovState = cellDrMapping[l_face].godunov;
      drKrnl.Qext = data.dofs;
      drKrnl._prefetch.I = faceNeighbors_prefetch[l_face];
      drKrnl.execute(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
    }
  }

  kernel::neighbour nKrnl = m_nKrnlPrototype;
  nKrnl.Qext = Qext;
  nKrnl.Q = data.dofs;
  nKrnl.Qane = data.dofsAne;
  nKrnl.w = data.neighboringIntegration.specific.w;
  
  nKrnl.execute();
}

void seissol::kernels::Neighbor::flopsNeighborsIntegral( const enum faceType  i_faceTypes[4],
                                                         const int            i_neighboringIndices[4][2],
                                                         CellDRMapping const (&cellDrMapping)[4],
                                                         unsigned int        &o_nonZeroFlops,
                                                         unsigned int        &o_hardwareFlops,
                                                         long long&           o_drNonZeroFlops,
                                                         long long&           o_drHardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;
  o_drNonZeroFlops = 0; o_drHardwareFlops = 0;
  
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != outflow && i_faceTypes[l_face] != dynamicRupture ) {
      // compute the neighboring elements flux matrix id.
      if( i_faceTypes[l_face] != freeSurface ) {
        assert(i_neighboringIndices[l_face][0] < 4 && i_neighboringIndices[l_face][1] < 3);
        
        o_nonZeroFlops  += seissol::kernel::neighbourFluxExt::nonZeroFlops(i_neighboringIndices[l_face][1], i_neighboringIndices[l_face][0], l_face);
        o_hardwareFlops += seissol::kernel::neighbourFluxExt::hardwareFlops(i_neighboringIndices[l_face][1], i_neighboringIndices[l_face][0], l_face);
      } else { // fall back to local matrices in case of free surface boundary conditions
        o_nonZeroFlops  += seissol::kernel::localFluxExt::nonZeroFlops(l_face);
        o_hardwareFlops += seissol::kernel::localFluxExt::hardwareFlops(l_face);
      }
    } else if (i_faceTypes[l_face] == dynamicRupture) {
      o_drNonZeroFlops += kernel::nodalFlux::nonZeroFlops(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
      o_drHardwareFlops += kernel::nodalFlux::hardwareFlops(cellDrMapping[l_face].side, cellDrMapping[l_face].faceRelation);
    }
  }

  o_nonZeroFlops += kernel::neighbour::NonZeroFlops;
  o_hardwareFlops += kernel::neighbour::HardwareFlops;
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
