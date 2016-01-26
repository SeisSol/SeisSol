/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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

#include "Neighbor.h"

#ifndef NDEBUG
#pragma message "compiling boundary kernel with assertions"
#endif

#include <generated_code/kernels.h>
#include <generated_code/flops.h>

#include <cassert>
#include <stdint.h>
#include <cstddef>

void seissol::kernels::Neighbor::computeNeighborsIntegral(  enum faceType const               i_faceTypes[4],
                                                            int const                         i_neighboringIndices[4][2],
                                                            GlobalData const*                 global,
                                                            NeighboringIntegrationData const* neighbor,
                                                            real*                             i_timeIntegrated[4],
                                                            real                              io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ] ) {
#ifndef NDEBUG
  // alignment of the flux matrices
  for( int l_matrix = 0; l_matrix < 52; l_matrix++ ) {
    assert( ((uintptr_t)global->fluxMatrices[l_matrix]) % ALIGNMENT == 0 );
  }

  // alignment of the time integrated dofs
  for( int l_neighbor = 0; l_neighbor < 4; l_neighbor++ ) {
    if( i_faceTypes[l_neighbor] != outflow && i_faceTypes[l_neighbor] != dynamicRupture ) { // no alignment for outflow and DR boundaries required
      assert( ((uintptr_t)i_timeIntegrated[l_neighbor]) % ALIGNMENT == 0 );
    }
  }
#endif

  // alignment of the degrees of freedom
  assert( ((uintptr_t)io_degreesOfFreedom) % ALIGNMENT == 0 );

  /*
   * compute neighboring cells contribution to the boundary integral.
   */
  // iterate over faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != outflow && i_faceTypes[l_face] != dynamicRupture ) {
      // compute the neighboring elements flux matrix id.
      if( i_faceTypes[l_face] != freeSurface ) {
        // derive memory and kernel index
        unsigned l_id = l_face*12                          // jump over index \f$i\f$
                      + i_neighboringIndices[l_face][0]*3  // jump over index \f$j\f$
                      + i_neighboringIndices[l_face][1];   // jump over index \f$h\f$
        
        // assert we have a valid index.
        assert( l_id < 48 );
        
        seissol::generatedKernels::neighboringFlux[l_id](
          neighbor->nAmNm1[l_face],
          global->fluxMatrices[4 + l_id],
          i_timeIntegrated[l_face],
          io_degreesOfFreedom
        );
      } else { // fall back to local matrices in case of free surface boundary conditions
        seissol::generatedKernels::localFlux[l_face](
          neighbor->nAmNm1[l_face],
          global->fluxMatrices[l_face],
          i_timeIntegrated[l_face],
          io_degreesOfFreedom
        );
      }
    }
  }
}

void seissol::kernels::Neighbor::flopsNeighborsIntegral( const enum faceType  i_faceTypes[4],
                                                         const int            i_neighboringIndices[4][2],
                                                         unsigned int        &o_nonZeroFlops,
                                                         unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;
  
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != outflow && i_faceTypes[l_face] != dynamicRupture ) {
      // compute the neighboring elements flux matrix id.
      if( i_faceTypes[l_face] != freeSurface ) {
        // derive memory and kernel index
        unsigned l_id = l_face*12                          // jump over index \f$i\f$
                      + i_neighboringIndices[l_face][0]*3  // jump over index \f$j\f$
                      + i_neighboringIndices[l_face][1];   // jump over index \f$h\f$
        
        // assert we have a valid index.
        assert( l_id < 48 );
        
        o_nonZeroFlops  += seissol::flops::neighboringFlux_nonZero[l_id];
        o_hardwareFlops += seissol::flops::neighboringFlux_hardware[l_id];
      } else { // fall back to local matrices in case of free surface boundary conditions
        o_nonZeroFlops  += seissol::flops::localFlux_nonZero[l_face];
        o_hardwareFlops += seissol::flops::localFlux_hardware[l_face];
      }
    }
  }
}

unsigned seissol::kernels::Neighbor::bytesNeighborsIntegral()
{
  unsigned reals = 0;

  // 4 * tElasticDOFS load, DOFs load, DOFs write
  reals += 4 * NUMBER_OF_ALIGNED_ELASTIC_DOFS + 2 * NUMBER_OF_ALIGNED_DOFS;
  // flux solvers load
  reals += 4 * seissol::model::AminusT::reals;
  
  return reals * sizeof(real);
}
