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
 * Volume kernel of SeisSol.
 **/

#ifndef VOLUME_H_
#define VOLUME_H_

#include "typedefs.hpp"
#include <cassert>
#include "common.hpp"

namespace seissol {
  namespace kernels {
    class Volume;
  }
}

/**
 * Volume kernel, which computes the volume integration.
 **/
class seissol::kernels::Volume {
  // explicit private for unit tests
  private:
    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global stiffness matrix (case a) or B the sparse star matrix (case b).
     * Each kernel can be dense or sparse, the kernels are ordered as follows:
     *    0:  \f$ K^\xi   \f$
     *    1:  \f$ K^\eta  \f$
     *    2:  \f$ K^\zeta \f$
     *    3:  \f$ A^* \vee B^* \vee C^* \f
     *
     * The matrix kernels might prefetch matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     * 
     * @param i_A left/stiffness matrix (case a) or time integrated DOFs matrix (case b).
     * @param i_B right/time integrated DOFs matrix (case a) or star matrix (case b).
     * @param io_C result matrix.
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[4])( const real *i_A,         const real *i_B,               real *io_C,
                                const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );

    /**
     * Number of non-zero floating point operations performed by each matrix kernel.
     **/
    unsigned int m_nonZeroFlops[4];

    /**
     * Number of floating point operations in hardware performed by each matrix kernels
     **/
    unsigned int m_hardwareFlops[4];

  public:
    /**
     * Constructor, which initializes the volume kernel.
     **/
    Volume();

    /**
     * Computes the volume integral from previously computed time integrated degrees of freedom.
     *
     * @param i_stiffnessMatrices stiffness matrices, 0: \f$ K^\xi \f$, 1: \f$ K^\eta\f$, 2: \f K^\zeta \f$.
     * @param i_timeIntegratedDegreesOfFreedom time integrated degrees of freedom.
     * @param i_starMatrices star matrices, 0: \f$ A^*_k \f$, 1: \f$ B^*_k \f$, 2: \f$ C^*_k \f$.
     * @param io_degreesOfFreedom degrees of freedom.
     **/
    void computeIntegral( real** i_stiffnessMatrices,
                          real*  i_timeIntegratedDegreesOfFreedom,
                          real   i_starMatrices[3][STAR_NNZ],
                          real*  io_degreesOfFreedom );

    /**
     * Derives the number of non-zero and hardware floating point operation in the volume integration.
     * @param o_nonZeroFlops number of performed non zero floating point operations.
     * @param o_hardwareFlops number of performed floating point operations in hardware.
     **/
    void flopsIntegral( unsigned int &o_nonZeroFlops,
                        unsigned int &o_hardwareFlops );
};

#endif

