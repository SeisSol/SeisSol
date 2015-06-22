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

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "typedefs.hpp"

namespace seissol {
  namespace kernels {
    class Boundary;
  }
}

/**
 * Boundary/Flux kernel, which computes the boundary integration.
 **/
class seissol::kernels::Boundary {
  // explicit private for unit tests
  private:
    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global flux matrix (case a) or B the flux solver (case b, 52).
     *   Each kernel can be dense or sparse.
     *   The kernels are ordered element local contribution \f$ F^{-, i} \f$ in front, which is 
     *   followed by the flux matrices for the neighbor element contribution \f$ F^{+, i, j, h} \f$
     *    0-51:  \f$ F^{-, 1} \vee \ldots \vee F^{-, 4} \vee F^+{+, 1, 1, 1} \vee \ldots \vee F^+{+, 4, 4, 3} \f$
     *    52:    \f$ N_{k,i} A_k^+ N_{k,i}^{-1}\f$ or \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$ without prefetches
     *    53:    \f$ N_{k,i} A_k^+ N_{k,i}^{-1}\f$ or \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$ with prefetches
     *
     * The matrix kernels might prefetch matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     *
     * @param i_A left/flux matrix (case a) or unknowns matrix (case b).
     * @param i_B right/unknowns matrix (case a) or flux solver (case b).
     * @param io_C result matrix.
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[54])( const real *i_A,         const real *i_B,               real *io_C,
                                 const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );

    /**
     * Number of non-zero floating point operations performed by each matrix kernel.
     **/
    unsigned int m_nonZeroFlops[54];

    /**
     * Number of floating point operations in hardware performed by each matrix kernels
     **/
    unsigned int m_hardwareFlops[54];

  public:
    /**
     * Constructor, which initializes the boundary kernel.
     **/
    Boundary();

    /**
     * Computes the cell's local contribution to the boundary integral.
     *
     * @param i_faceTypes types of the faces: regular, freeSurface, dynamicRupture or periodic
     * @param i_fluxMatrices 52 flux matrices in the following order:
     *        0:  \f$ F^{-, 1} \f$
     *        1 : \f$ F^{-, 2} \f$
     *        2:  \f$ F^{-, 3} \f$
     *        3 : \f$ F^{-, 4} \f$
     *        4:  \f$ F^+{+, 1, 1, 1} \f$
     *        5:  \f$ F^+{+, 1, 1, 2} \f$
     *        6:  \f$ F^+{+, 1, 1, 3} \f$
     *        7:  \f$ F^+{+, 1, 2, 1} \f$
     *        8:  \f$ F^+{+, 1, 1, 2} \f$
     *        9:  \f$ F^+{+, 1, 1, 3} \f$
     *        [..]
     *        51: \f$ F^+{+, 4, 4, 3} \f$.
     * @param i_timeIntegrated time integrated degrees of freedoms of the cell \f$i \f$.
     * @param i_fluxSolvers matrices \f$N_{k,i} A_k^+ N_{k,i}^{-1}\f$,
     *        where \f$N_{k,i}^{-1}\f$ and \f$N_{k,i}\f$ account for the forth and back transformation
     *        to reference space relative to the outerpointing normal of face \f$i\f$ in element
     *        \f$k\f$. \f$A_k^+\f$ is the flux contribution of cell \f$k\f$, thus assembled using
     *        the positive eigenvalues of element \f$k\f$.
     * @param io_degreesOfFreedom DOFs, which will be updated by the boundary integral.
     **/
    void computeLocalIntegral( const enum faceType i_faceTypes[4],
                                     real         *i_fluxMatrices[52],
                                     real          i_timeIntegrated[    NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ],
                                     real          i_fluxSolvers[4][    NUMBER_OF_QUANTITIES             *NUMBER_OF_QUANTITIES ],
                                     real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ] );

    /**
     * Derives the number of sparse and hardware floating point operations performed in the local boundary integral.
     *
     * @param i_faceTypes types of the faces: regular, freeSurface, dynamicRupture or periodic.
     * @param o_nonZeroFlops will be set to the number of non-zero flops performed by the element-local contribution to the boundary integration.
     * @param o_hardwareFlops will be set to the number of flops in hardware performed by the element-local contribution to the boundary integration.
     **/
    void flopsLocalIntegral( const enum faceType  i_faceTypes[4],
                             unsigned int        &o_nonZeroFlops,
                             unsigned int        &o_hardwareFlops );

    /**
     * Computes the neighboring cells contribution to the boundary integral for a single cell.
     *
     * @param i_faceTypes types of the faces: regular, freeSurface, dynamicRupture or periodic
     * @param i_neighboringIndices indices \f$j(i) \in \{0,1,2,3\}\f$ and \f$h(i) \in \{0,1,2\}\f$,
     *        which depend on the face combinations of the current elements faces \f$ i \in
     *        \{0,1,2,3\}\f$ and the neighboring faces \f$ j(i) \f$ and vertex orientation of the
     *        element and neighboring face \f$ h(i) \f$.
     *        A two dimensional array is expected, where the first index i_neighboringIndices[i][]
     *        selects the face of the current element and the second index i_neighboringIndices[][*]
     *        gives \f$ j(i)\f$: i_neighboringIndices[][0] or \f$ h(i) \f$:  i_neighboringIndices[][1].
     * @param i_fluxMatrices 52 flux matrices in the following order:
     *        0:  \f$ F^{-, 1} \f$
     *        1 : \f$ F^{-, 2} \f$
     *        2:  \f$ F^{-, 3} \f$
     *        3 : \f$ F^{-, 4} \f$
     *        4:  \f$ F^+{+, 1, 1, 1} \f$
     *        5:  \f$ F^+{+, 1, 1, 2} \f$
     *        6:  \f$ F^+{+, 1, 1, 3} \f$
     *        7:  \f$ F^+{+, 1, 2, 1} \f$
     *        8:  \f$ F^+{+, 1, 1, 2} \f$
     *        9:  \f$ F^+{+, 1, 1, 3} \f$
     *        [..]
     *        51: \f$ F^+{+, 4, 4, 3} \f$.
     * @param i_timeIntegrated time integrated degrees of freedoms of the neighboring cells \f$k(i) \f$.
     * @param i_fluxSolvers flux solver matrices \f$N_{k,i} A_{k(i)}^- N_{k,i}^{-1}\f$,
     *        where \f$N_{k,i}^{-1}\f$ and \f$N_{k,i}\f$ account for the forth and back
     *        transformation to reference space relative to the outerpointing normal of face $i$ in
     *        element \f$k\f$. \f$A_{k(i)}^-\f$ is the flux contribuation of the \f$i\f$-th face
     *        neighbor of cell \f$k\f$, thus assembled using negative eigenvalues of element \f$k(i)\f$.
     * @param io_degreesOfFreedom DOFs, which will be updated by the boundary integral.
     * @param i_faceNeighbors_prefetch time integrated degrees of freedoms/time derivates of the neighboring cells 
     *        which should be fetched into L2 cache in this call. There is no garantuee that these prefetchings
     *        are acutally going to be executed. It depends on the employed GEMM implementation.
     * @param i_fluxMatricies_prefetch pointer to a flux matrix which should be fetched into L2 cache in this call. 
     *        There is no garantuee that these prefetchings are acutally going to be executed. 
     *        It depends on the employed GEMM implementation.
     **/
    void computeNeighborsIntegral( const enum faceType i_faceTypes[4],
                                   const int           i_neighboringIndices[4][2],
                                         real         *i_fluxMatrices[52],
                                         real         *i_timeIntegrated[4],
                                         real          i_fluxSolvers[4][    NUMBER_OF_QUANTITIES             *NUMBER_OF_QUANTITIES ],
#ifdef ENABLE_MATRIX_PREFETCH
                                         real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ],
                                         real         *i_faceNeighbors_prefetch[4],
                                         real         *i_fluxMatricies_prefetch[4] );
#else
                                         real          io_degreesOfFreedom[ NUMBER_OF_ALIGNED_BASIS_FUNCTIONS*NUMBER_OF_QUANTITIES ] );
#endif

    /**
     * Derives the number of sparse and hardware floating point opreations performed in the neighboring intergal.
     * @param i_faceTypes types of the faces: regular, freeSurface, dynamicRupture or periodic
     * @param i_neighboringIndices indices \f$j(i) \in \{0,1,2,3\}\f$ and \f$h(i) \in \{0,1,2\}\f$,
     *        which depend on the face combinations of the current elements faces \f$ i \in
     *        \{0,1,2,3\}\f$ and the neighboring faces \f$ j(i) \f$ and vertex orientation of the
     *        element and neighboring face \f$ h(i) \f$.
     *        A two dimensional array is expected, where the first index i_neighboringIndices[i][]
     *        selects the face of the current element and the second index i_neighboringIndices[][*]
     * @param o_nonZeroFlops will be set to the number of non-zero flops performed by the contribution of the neighboring cells to the boundary integration.
     * @param o_hardwareFlops will be set to the number of flops in hardware performed by the contribution of the neighboring cells to the boundary integration.
     **/
    void flopsNeighborsIntegral( const enum faceType  i_faceTypes[4],
                                 const int            i_neighboringIndices[4][2],
                                 unsigned int        &o_nonZeroFlops,
                                 unsigned int        &o_hardwareFlops );
};

#endif
