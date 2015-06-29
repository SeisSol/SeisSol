/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
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
 * Common kernel-level functions
 **/

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <algorithm>
#include <Initializer/typedefs.hpp>

namespace seissol {
  namespace kernels {
    /**
     * Gets the number of basis functions for the given convergence order.
     *
     * @param i_convergenceOrder convergence order.
     * @return number of basis funcitons.
     **/
    static unsigned int getNumberOfBasisFunctions( unsigned int i_convergenceOrder = CONVERGENCE_ORDER ) {
      return i_convergenceOrder*(i_convergenceOrder+1)*(i_convergenceOrder+2)/6;
    }

    /**
     * Gets the number of aligned reals.
     *
     * @param i_alignment alignment in bytes.
     * @return aligned number of reals.
     **/
    static unsigned int getNumberOfAlignedReals( unsigned int i_numberOfReals,
                                                 unsigned int i_alignment = ALIGNMENT ) {

      unsigned int l_nonZeroBytes = i_numberOfReals * sizeof(real);

      // compute corresponding aligned # of basis functions
      unsigned int l_alignedBytes = l_nonZeroBytes + (i_alignment-(l_nonZeroBytes % i_alignment))%i_alignment;
      assert( l_alignedBytes % sizeof(real) == 0 );

      return l_alignedBytes / sizeof( real );
    }

    /**
     * Get the # of basis functions aligned to the given boundaries.
     *
     * @param i_convergenceOrder convergence order.
     * @param i_alignment alignment in bytes.
     * @return aligned number of basis functions.
     **/
    static unsigned int getNumberOfAlignedBasisFunctions( unsigned int i_convergenceOrder = CONVERGENCE_ORDER,
                                                          unsigned int i_alignment        = ALIGNMENT ) {
      unsigned int l_numberOfBasisFunctions = getNumberOfBasisFunctions( i_convergenceOrder);

      return getNumberOfAlignedReals( l_numberOfBasisFunctions );
    }

    /**
     * Converts memory aligned degrees of freedom (with zero padding) to unaligned (compressed, without zero padding) storage.
     *
     * @param i_alignedDofs aligned degrees of freedom (zero padding in the basis functions / columns).
     * @param o_unalignedDofs unaligned degrees of freedom.
     **/
    template<typename real_from, typename real_to>
    static void convertAlignedDofs( const real_from i_alignedDofs[   NUMBER_OF_ALIGNED_DOFS],
                                          real_to   o_unalignedDofs[ NUMBER_OF_DOFS] ) {
      for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
        for( unsigned int l_basisFunction = 0; l_basisFunction < NUMBER_OF_BASIS_FUNCTIONS; l_basisFunction++ ) {
          o_unalignedDofs[l_quantity*NUMBER_OF_BASIS_FUNCTIONS + l_basisFunction] = i_alignedDofs[l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + l_basisFunction];
        }
      }
    }

    /**
     * Converts unaligned degrees of freedom (compressed, without zero padding) to aligned storage (with zero padding).
     *   Remark: It is assumed that the overhead of the alignment is initialized with zero.
     *
     * @param i_unalignedDofs unaligned degrees of freedom.
     * @param o_alignedDofs aligned degrees of freedom (zero paddin in the basis functions / columns).
     **/
    template<typename real_from, typename real_to>
    static void convertUnalignedDofs( const real_from i_unalignedDofs[ NUMBER_OF_DOFS ],
                                            real_to   o_alignedDofs[   NUMBER_OF_ALIGNED_DOFS ] ) {
      for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
        for( unsigned int l_basisFunction = 0; l_basisFunction < NUMBER_OF_BASIS_FUNCTIONS; l_basisFunction++ ) {
          o_alignedDofs[l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + l_basisFunction] = i_unalignedDofs[l_quantity*NUMBER_OF_BASIS_FUNCTIONS + l_basisFunction];
        }
      }
    }

    /**
     * Adds the unaligned update to degrees of freedom with aligned storage.
     *
     * @param i_unalignedUpdate unaligned update.
     * @param o_alignedDofs aligned degrees of freedom.
     **/
    template<typename real_from, typename real_to>
    static void addToAlignedDofs( const real_from i_unalignedUpdate[ NUMBER_OF_DOFS ],
                                        real_to   o_alignedDofs[     NUMBER_OF_ALIGNED_DOFS ] ) {
      for( unsigned int l_quantity = 0; l_quantity < NUMBER_OF_QUANTITIES; l_quantity++ ) {
        for( unsigned int l_basisFunction = 0; l_basisFunction < NUMBER_OF_BASIS_FUNCTIONS; l_basisFunction++ ) {
          o_alignedDofs[l_quantity*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + l_basisFunction] += i_unalignedUpdate[l_quantity*NUMBER_OF_BASIS_FUNCTIONS + l_basisFunction];
        }
      }
    }

    /*
     * Copies a submatrix of A (sizes of B) to B.
     * If B doesn't fit in A zeros are set.
     *
     * @param i_A values of matrix A.
     * @param i_aNumberOfRows number of rows of A.
     * @param i_aNumberOfColumns number of columns of A.
     * @param i_aLeadingDimension leading dimension of A.
     * @param o_B values of matrix B, which will be set.
     * @param i_bNumberOfRows number of rows of matrix B.
     * @param i_bNumberOfColumns number of columns of matrix B.
     * @param i_bLeadingDimension leading dimension of B.
     */
    template<typename real_from, typename real_to>
    static void copySubMatrix( const real_from    *i_A,
                               const unsigned int  i_aNumberOfRows,
                               const unsigned int  i_aNumberOfColumns,
                               const unsigned int  i_aLeadingDimension,
                                     real_to      *o_B,
                               const unsigned int  i_bNumberOfRows,
                               const unsigned int  i_bNumberOfColumns,
                               const unsigned int  i_bLeadingDimension ) {
      // set matrix B to zero
      for( unsigned int l_index = 0; l_index < i_bLeadingDimension*i_bNumberOfColumns; l_index++ ) {
        o_B[l_index] = (real) 0;
      }

      // copy the entries
      for( unsigned int l_column = 0; l_column < std::min( i_aNumberOfColumns, i_bNumberOfColumns ); l_column++ ) {
        for( unsigned int l_row = 0; l_row < std::min( i_aNumberOfRows, i_bNumberOfRows ); l_row++ ) {
          unsigned int l_aIndex = l_column * i_aLeadingDimension + l_row;
          unsigned int l_bIndex = l_column * i_bLeadingDimension + l_row;

          o_B[l_bIndex] = i_A[l_aIndex];
        }
      }
    }

    /**
     * Convert compressed and memory aligned time derivatives to a full (including zeros) unaligned format.
     *
     * @param i_compressedDerivatives derivatives in compressed, aligned format.
     * @param o_fullDerivatives derivatives in full, unaligned format.
     **/
    template<typename real_from, typename real_to>
    static void convertAlignedCompressedTimeDerivatives( const real_from *i_compressedDerivatives,
                                                               real_to    o_fullDerivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] ) {
      unsigned int l_firstEntry = 0;

      for( unsigned int l_order = 0; l_order < CONVERGENCE_ORDER; l_order++ ) {
        copySubMatrix( &i_compressedDerivatives[l_firstEntry],
                        getNumberOfBasisFunctions( CONVERGENCE_ORDER-l_order ),
                        NUMBER_OF_QUANTITIES,
                        getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order ),
                        o_fullDerivatives[l_order],
                        NUMBER_OF_BASIS_FUNCTIONS,
                        NUMBER_OF_QUANTITIES,
                        NUMBER_OF_BASIS_FUNCTIONS );

        l_firstEntry += getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-l_order ) * NUMBER_OF_QUANTITIES;
      }
    }

    /**
     * Converts a full star matrix (including zeros) to a compressed star matrix containing the 24 possible non-zeros.
     *
     * @param i_fullStarMatrix star matrix in full format.
     * @param o_compressedStarMatrix compressed star matrix.
     **/
    template<typename real_from, typename real_to>
    static void convertStarMatrix( const real_from *i_fullStarMatrix,
                                         real_to   *o_compressedStarMatrix ) {
      o_compressedStarMatrix[ 0] = i_fullStarMatrix[0*9 + 6];
      o_compressedStarMatrix[ 1] = i_fullStarMatrix[0*9 + 7];
      o_compressedStarMatrix[ 2] = i_fullStarMatrix[0*9 + 8];

      o_compressedStarMatrix[ 3] = i_fullStarMatrix[1*9 + 6];
      o_compressedStarMatrix[ 4] = i_fullStarMatrix[1*9 + 7];
      o_compressedStarMatrix[ 5] = i_fullStarMatrix[1*9 + 8];

      o_compressedStarMatrix[ 6] = i_fullStarMatrix[2*9 + 6];
      o_compressedStarMatrix[ 7] = i_fullStarMatrix[2*9 + 7];
      o_compressedStarMatrix[ 8] = i_fullStarMatrix[2*9 + 8];

      o_compressedStarMatrix[ 9] = i_fullStarMatrix[3*9 + 6];
      o_compressedStarMatrix[10] = i_fullStarMatrix[3*9 + 7];

      o_compressedStarMatrix[11] = i_fullStarMatrix[4*9 + 7];
      o_compressedStarMatrix[12] = i_fullStarMatrix[4*9 + 8];

      o_compressedStarMatrix[13] = i_fullStarMatrix[5*9 + 6];
      o_compressedStarMatrix[14] = i_fullStarMatrix[5*9 + 8];

      o_compressedStarMatrix[15] = i_fullStarMatrix[6*9 + 0];
      o_compressedStarMatrix[16] = i_fullStarMatrix[6*9 + 3];
      o_compressedStarMatrix[17] = i_fullStarMatrix[6*9 + 5];

      o_compressedStarMatrix[18] = i_fullStarMatrix[7*9 + 1];
      o_compressedStarMatrix[19] = i_fullStarMatrix[7*9 + 3];
      o_compressedStarMatrix[20] = i_fullStarMatrix[7*9 + 4];

      o_compressedStarMatrix[21] = i_fullStarMatrix[8*9 + 2];
      o_compressedStarMatrix[22] = i_fullStarMatrix[8*9 + 4];
      o_compressedStarMatrix[23] = i_fullStarMatrix[8*9 + 5];
    }

  }
}

#endif
