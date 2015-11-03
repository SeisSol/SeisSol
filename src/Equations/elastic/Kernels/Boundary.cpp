/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
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

#include "Boundary.h"

#include <matrix_kernels/sparse.h>
#include <matrix_kernels/dense.h>

#ifdef ENABLE_STREAM_MATRIX_PREFETCH
#include <immintrin.h>
#endif

#ifndef NDEBUG
#pragma message "compiling boundary kernel with assertions"
#endif

#include <cassert>
#include <stdint.h>
#include <cstddef>

seissol::kernels::Boundary::Boundary() {
  // intialize the function pointers to the matrix kernels
#define BOUNDARY_KERNEL
#include <initialization/bind.h>
#undef BOUNDARY_KERNEL
}

void seissol::kernels::Boundary::computeLocalIntegral( const enum faceType i_faceTypes[4],
                                                             real         *i_fluxMatrices[52],
                                                             real          i_timeIntegrated[    NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ],
                                                             real          i_fluxSolvers[4][    NUMBER_OF_QUANTITIES             *NUMBER_OF_QUANTITIES ],
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
                                                             real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ],
                                                             real          i_timeIntegrated_prefetch[    NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ],
                                                             real          io_degreesOfFreedom_prefetch[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ] ) {
#else
                                                             real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ] ) {
#endif
  /*
   * assert valid input
   */
#ifndef NDEBUG
  // alignment of the flux matrices
  for( int l_matrix = 0; l_matrix < 52; l_matrix++ ) {
    assert( ((uintptr_t)i_fluxMatrices[l_matrix]) % ALIGNMENT == 0 );
  }
#endif
  // alignment of the time integrated dofs
  assert( ((uintptr_t)i_timeIntegrated) % ALIGNMENT == 0 );

  // alignment of the degrees of freedom
  assert( ((uintptr_t)io_degreesOfFreedom) % ALIGNMENT == 0 );

  /*
   * compute cell local contribution of the boundary integral.
   */
  // temporary product (we have to multiply a matrix from the left and the right)
  real l_temporaryResult[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));

  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != dynamicRupture ) {
      // compute neighboring elements contribution
      m_matrixKernels[l_face]( i_fluxMatrices[l_face], i_timeIntegrated,      l_temporaryResult,
                               NULL,                   NULL,                  NULL                 ); // These will be be ignored
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
      if (l_face == 1) {
        m_matrixKernels[54](     l_temporaryResult,      i_fluxSolvers[l_face],       io_degreesOfFreedom,
                                 NULL,                   i_timeIntegrated_prefetch,   NULL                 );
      } else if (l_face == 2) {
        m_matrixKernels[54](     l_temporaryResult,      i_fluxSolvers[l_face],              io_degreesOfFreedom,
                                 NULL,                   io_degreesOfFreedom_prefetch,       NULL                 );
      } else {
        m_matrixKernels[52](     l_temporaryResult,      i_fluxSolvers[l_face], io_degreesOfFreedom,
                                 NULL,                   NULL,                  NULL                 ); // These will be be ignored
      }
#else
      m_matrixKernels[52](     l_temporaryResult,      i_fluxSolvers[l_face], io_degreesOfFreedom,
                               NULL,                   NULL,                  NULL                 ); // These will be be ignored
#endif
    } else {
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
      if (l_face == 0) {
        for (unsigned int l_p = 0; l_p < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES; l_p += 8 ){
          _mm_prefetch( (const char*) i_timeIntegrated_prefetch+l_p, _MM_HINT_T1 );
        }
      } else if (l_face == 1) {
        for (unsigned int l_p = 0; l_p < NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES; l_p += 8 ){
          _mm_prefetch( (const char*) io_degreesOfFreedom_prefetch+l_p, _MM_HINT_T1 );
        }
      }
#endif
    }
  }
}

void seissol::kernels::Boundary::flopsLocalIntegral( const enum faceType  i_faceTypes[4],
                                                     unsigned int        &o_nonZeroFlops,
                                                     unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;

  // iterate over faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != dynamicRupture ) {
      // compute number of flops
      o_nonZeroFlops  += m_nonZeroFlops[ l_face];
      o_hardwareFlops += m_hardwareFlops[l_face];

      o_nonZeroFlops  += m_nonZeroFlops[ 52];
      o_hardwareFlops += m_hardwareFlops[52];
    }
  }
}

void seissol::kernels::Boundary::computeNeighborsIntegral( const enum faceType i_faceTypes[4],
                                                           const int           i_neighboringIndices[4][2],
                                                                 real         *i_fluxMatrices[52],
                                                                 real         *i_timeIntegrated[4],
                                                                 real          i_fluxSolvers[4][    NUMBER_OF_QUANTITIES             *NUMBER_OF_QUANTITIES ],
#ifdef ENABLE_MATRIX_PREFETCH
                                                                 real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ],
                                                                 real         *i_faceNeighbors_prefetch[4],
                                                                 real         *i_fluxMatricies_prefetch[4] ) {
#else
                                                                 real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ] ) {
#endif

  /*
   * assert valid input
   */
#ifndef NDEBUG
  // alignment of the flux matrices
  for( int l_matrix = 0; l_matrix < 52; l_matrix++ ) {
    assert( ((uintptr_t)i_fluxMatrices[l_matrix]) % ALIGNMENT == 0 );
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
  // temporary product (we have to multiply a matrix from the left and the right)
  real l_temporaryResult[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(PAGESIZE_STACK)));

#ifdef ENABLE_STREAM_MATRIX_PREFETCH
  // perfetch the first cacheline per column for the DOFs
  for ( unsigned int l_var = 0; l_var < NUMBER_OF_QUANTITIES; l_var++ ) {
    _mm_prefetch( (const char*) io_degreesOfFreedom+(l_var*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS), _MM_HINT_T1 );
  }
#endif

  // iterate over faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != outflow && i_faceTypes[l_face] != dynamicRupture ) {
      // id of the flux matrix (0-3: element local, 4-51: neighboring element)
      unsigned int l_id;

      // compute the neighboring elements flux matrix id.
      if( i_faceTypes[l_face] != freeSurface ) {
        // derive memory and kernel index
        l_id = 4                                   // jump over flux matrices \f$ F^{-, i} \f$
              + l_face*12                          // jump over index \f$i\f$
              + i_neighboringIndices[l_face][0]*3  // jump over index \f$j\f$
              + i_neighboringIndices[l_face][1];   // jump over index \f$h\f$

        // assert we have a neighboring index in the case of non-absorbing boundary conditions.
        assert( l_id >= 4 || i_faceTypes[l_face] == outflow );
      }
      else { // fall back to local matrices in case of free surface boundary conditions
        l_id = l_face;
      }

      // assert we have a valid index.
      assert( l_id < 52 );

#ifdef ENABLE_MATRIX_PREFETCH
      real* l_faceNeigh_prefetch = NULL;
      real* l_flux_prefetch = NULL;
      if (i_faceNeighbors_prefetch[l_face] != NULL) {
        l_faceNeigh_prefetch = i_faceNeighbors_prefetch[l_face];
      } else {
        l_faceNeigh_prefetch = io_degreesOfFreedom;
      }
      if (i_fluxMatricies_prefetch[l_face] != NULL) {
        l_flux_prefetch = i_fluxMatricies_prefetch[l_face];
      } else {
        l_flux_prefetch = i_fluxMatrices[l_id];
      }

      // compute neighboring elements contribution
      m_matrixKernels[l_id]( i_fluxMatrices[l_id], i_timeIntegrated[l_face], l_temporaryResult,
                             l_flux_prefetch,      l_faceNeigh_prefetch,     NULL                 ); 

      m_matrixKernels[53](   l_temporaryResult,    i_fluxSolvers[l_face],    io_degreesOfFreedom,
                             l_flux_prefetch,      l_faceNeigh_prefetch,     NULL                 ); 
#else
      // compute neighboring elements contribution
      m_matrixKernels[l_id]( i_fluxMatrices[l_id], i_timeIntegrated[l_face], l_temporaryResult,
                             NULL,                 NULL,                     NULL                 );  // These will be be ignored

      m_matrixKernels[53](   l_temporaryResult,    i_fluxSolvers[l_face],    io_degreesOfFreedom,
                             NULL,                 NULL,                     NULL                 );  // These will be be ignored
#endif
    }
  }
}

void seissol::kernels::Boundary::flopsNeighborsIntegral( const enum faceType  i_faceTypes[4],
                                                         const int            i_neighboringIndices[4][2],
                                                         unsigned int        &o_nonZeroFlops,
                                                         unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops = 0;

  // iterate over faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // no neighboring cell contribution in the case of absorbing and dynamic rupture boundary conditions
    if( i_faceTypes[l_face] != outflow && i_faceTypes[l_face] != dynamicRupture ) {
      // id of the flux matrix (0-3: element local, 4-51: neighboring element)
      unsigned int l_id;

      // compute the neighboring elements flux matrix id.
      if( i_faceTypes[l_face] != freeSurface ) {
        // derive memory and kernel index
        l_id = 4                                   // jump over flux matrices \f$ F^{-, i} \f$
              + l_face*12                          // jump over index \f$i\f$
              + i_neighboringIndices[l_face][0]*3  // jump over index \f$j\f$
              + i_neighboringIndices[l_face][1];   // jump over index \f$h\f$

        // assert we have a neighboring index in the case of non-absorbing boundary conditions.
        assert( l_id >= 4 || i_faceTypes[l_face] == outflow );
      }
      else { // fall back to local matrices in case of free surface boundary conditions
        l_id = l_face;
      }

      // assert we have a valid index.
      assert( l_id < 52 );

      // compute number of flops
      o_nonZeroFlops  += m_nonZeroFlops[ l_id];
      o_hardwareFlops += m_hardwareFlops[l_id];

      o_nonZeroFlops  += m_nonZeroFlops[ 53];
      o_hardwareFlops += m_hardwareFlops[53];
    }
  }
}
