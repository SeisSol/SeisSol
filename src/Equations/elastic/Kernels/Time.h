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
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#ifndef TIME_H_
#define TIME_H_

#include <cassert>
#include <limits>
#include <Initializer/typedefs.hpp>
#include <Kernels/common.hpp>

namespace seissol {
  namespace kernels {
    class Time;
  }
}

/**
 * Time kernel, which computes the time derivatives, intgrals and extrapolation for an element using the ADER scheme.
 *
 * All functions operate on a compressed memory format, which stores only the (memory aligned) non-zeros of the time derivatives.
 *
 * Numerical motivation: The size of the non-zero time differentiated unknowns reduces in every step.
 *   Sketch for polynomial degree 2:
 * <pre>
 *   First step:
 *       \f$(k^\Xi)^T\f$     * \f$\frac{\partial^0}{\partial t^0} Q_k\f$ * \f$A^*_k\f$ + [...]
 *    _                   _     _                 _     _                 _             _                 _
 *   | - * - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * - - |           | * * * * * * * * * | <----/ remaining non-zero block
 *   | - - - - * - - - - - |   | * * * * * * * * * |   | - - - - - - - * - |           | * * * * * * * * * | <---/  size:
 *   | - - - - - * - - - - |   | * * * * * * * * * |   | - - - - - - - - * |           | * * * * * * * * * | <--/     #basis functions of polynomial degree 1
 *   | - - - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * * - |           | * * * * * * * * * | <-/
 *   | - - - - - - - - - - | * | * * * * * * * * * | * | - - - - - - - * * | + [...] = | - - - - - - - - - |    Remark: The star matrix multiplication is one
 *   | - - - - - - - - - - |   | * * * * * * * * * |   | - - - - - - * - * |           | - - - - - - - - - |            step ahead because of associativity
 *   | - - - - - - - - - - |   | * * * * * * * * * |   | * * * * - * - - - |           | - - - - - - - - - |            of the matrix multiplication.
 *   | - - - - - - - - - - |   | * * * * * * * * * |   | * * * * * - - - - |           | - - - - - - - - - |
 *   | - - - - - - - - - - |   | * * * * * * * * * |   |_* * * - * * - - -_|           |_- - - - - - - - -_|
 *   |_- - - - - - - - - -_|   |_* * * * * * * * *_|
 *
 *   Second step:
 *        \f$(k^\Xi)^T\f$     * \f$\frac{\partial^1}{\partial t^1} Q_k\f$ * \f$A^*_k\f$ + [...] 
 *    _                   _     _                 _     _                 _
 *   | - * - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * - - |
 *   | - - - - * - - - - - |   | * * * * * * * * * |   | - - - - - - - * - |
 *   | - - - - - * - - - - |   | * * * * * * * * * |   | - - - - - - - - * |
 *   | - - - - - * - * - - |   | * * * * * * * * * |   | - - - - - - * * - |
 *   | - - - - - - - - - - | * | - - - - - - - - - | * | - - - - - - - * * | + [...] = [...]
 *   | - - - - - - - - - - |   | - - - - - - - - - |   | - - - - - - * - * |
 *   | - - - - - - - - - - |   | - - - - - - - - - |   | * * * * - * - - - |
 *   | - - - - - - - - - - |   | - - - - - - - - - |   | * * * * * - - - - |
 *   | - - - - - - - - - - |   | - - - - - - - - - |   |_* * * - * * - - -_|
 *   |_- - - - - - - - - -_|   |_- - - - - - - - -_|
 *             ___________
 *            /     |     \
 *                  |
 *         columns, which hit zeros in the second step
 *      -> the product has only one non-zero row left
 * </pre>
 **/
class seissol::kernels::Time {
  // explicit private for unit tests
  private:
    //! aligned number of basis functions in decreasing order.
    unsigned int m_numberOfAlignedBasisFunctions[CONVERGENCE_ORDER];

    /*
     *! Offsets of the derivatives.
     *
     * * Offset counting starts at the zeroth derivative with o_derivativesOffset[0]=0; increasing derivatives follow:
     *   1st derivative: o_derivativesOffset[1]
     *   2nd derivative: o_derivativesOffset[2]
     *   ...
     * * Offset are always counted from positition zero; for example the sixth derivative will include all jumps over prior derivatives 0 to 5.
     */
    unsigned int m_derivativesOffsets[CONVERGENCE_ORDER];

    /**
     * Collection of matrix kernels, which perform the matrix product \f$ C += A.B\f$,
     * where \f$ A \f$ is a global transposed stiffness matrix (case a) or B a star matrix (case b).
     * Each matrix kernel can be dense or sparse.
     * Each kernel has hardcoded BLAS-specifiers (M, N, K, ld(A), ld(B), ld(C), beta) exploiting the recursive structure.
     * The kernels are ordered as follows:
     *    0-2:       1st derivative \f$ M^{-1} ( K^\xi )^T \vee M^{-1} ( K^\eta )^T \vee M^{-1} ( K^\zeta )^T \f$
     *    3:         1st derivative \f$ A^* \vee B^* \vee C^* \f
     *    4-6:       2nd derivative \f$ M^{-1} ( K^\xi )^T \vee M^{-1} ( K^\eta )^T \vee M^{-1} ( K^\zeta )^T \f$
     *    7:         2nd derivative \f$ A^* \vee B^* \vee C^* \f
     *    ...
     *    4*(O-2):   O-1th derivative
     *    4*(O-2)+1: O-1th derivative
     *
     * Remark: The mass matrix \f$ M \f$ is diagonal.
     * 
     * The matrix kernels might prefetch matrices of the next matrix multiplication triple \f$ A =+ B.C \f$,
     * thus loading upcoming matrices into lower level memory while the FPUs are busy.
     * In the case of the time integrator this means prefetching the transposed stiffness matrices (multiplied by the inverse mass matrices)
     * or star matrices of the upcoming operations in the recursive computation of the time derivatives.
     * The last operation,
     * \f[
     *   \left( \frac{\partial^{j_{\max-1}}}{\partial t^{j_{\max-1}} Q_k \right) C_k^*
     * \f],
     * already prefetches the stiffness, unknowns and time integrated unknowns matrix of the upcoming time integration of the next element.
     *
     * @param i_A left/transposed stiffness matrix (case a) or derivatives matrix (case b).
     * @param i_B right/derivatives matrix (case a) or star matrix (case b).
     * @param io_C resulting matrix.
     * @param i_APrefetch left matrix \f$ A \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_BPrefetch right matrix \f$ B \f$ of the next matrix triple \f$ (A, B, C) \f$.
     * @param i_CPrefetch result matrix \f$ C \f$ of the next matrix triple \f$ (A, B, C) \f$.
     **/  
    void (*m_matrixKernels[(CONVERGENCE_ORDER-1)*4])( const real *i_A,         const real *i_B,               real *io_C,
                                                      const real *i_APrefetch, const real *i_BPrefetch, const real *i_CPrefetch );

    /**
     * Number of non-zero floating point operations performed by each matrix kernel.
     **/
    unsigned int m_nonZeroFlops[(CONVERGENCE_ORDER-1)*4];

    /**
     * Number of floating point operations in hardware performed by each matrix kernels
     **/
    unsigned int m_hardwareFlops[(CONVERGENCE_ORDER-1)*4];

    /**
     * Stream-Copy (no-read-before-write) the degress of freedom into the first position of the time derivatives buffer
     *
     * @param i_degreesOfFreedom of the current time step \f$ t^\text{cell} \f$
     * @param o_derivativesBuffer time derivatives of the degrees of freedom in compressed format, this needs to be start address
     */
    inline void streamstoreFirstDerivative( const real*   i_degreesOfFreedom,
                                                  real*   o_derivativesBuffer );

    /**
     * Compute the time integraion of derivative i_derivative and add it to the timeIntegrated buffer o_timeIntegrated, optionally the 
     * the derivatives are stream-copied into o_timeDerivatives
     *
     * @param i_derivativeBuffer buffer containing the derivatives in compressed format
     * @param i_scalar the scalar factor of the time integration for the derivative which is currently being processed
     * @param i_derivative the current derivative
     * @param o_timeIntegrated the buffer into which the time integration is accumulated to
     * @param o_timeDerivatives (optional) time derivatives of the degrees of freedom in compressed format. If NULL only time integrated DOFs are returned.
     */
    inline void integrateInTime( const real*        i_derivativeBuffer,
                                       real         i_scalar,
                                       unsigned int i_derivative,
                                       real*        o_timeIntegrated,
                                       real*        o_timeDerivatives );

    /**
     * Initialize the timeIntegrated and derivativesBuffer before computing the time integration
     *
     * @param i_scalar the scalar factor of the time integration for the derivative which is currently being processed
     * @param i_degreesOfFreedom of the current time step \f$ t^\text{cell} \f$
     * @param o_timeIntegrated the buffer into which the time integration is accumulated to
     * @param o_derivativesBuffer time derivatives of the degrees of freedom in compressed format, this needs to be start address
     */
    inline void initialize( const real         i_scalar,
                            const real*        i_degreesOfFreedom,
                                  real*        o_timeIntegrated,
                                  real*        o_derivativesBuffer );

  public:
    /**
     * Constructor, which initializes the time kernel.
     **/
    Time();

    /**
     * Computes the ADER procedure.
     *   This is a hybrid call: output are the time integrated DOFs and (optional) the time derivatives.
     *   Storage format of the derivatives is compressed (storing only aligned non-zeros).
     *   Use computeTimeIntegral to compute time integrated degrees of freedom from the derivatives.
     *   Use computeTimeEvaluation to evaluate the time prediction of the degrees of freedom based on the derivatives.
     *
     * @param i_timeStepWidth time step width for the integration in time.
     * @param i_stiffnessMatrices negative transposed stiffness matrices (multiplied by inverse mass matrix), 0: \f$ -M^{-1} ( K^\xi )^T \f$ 1:\f$ -M^{-1} ( K^\eta )^T \f$ 2: \f$ -M^{-1} ( K^\zeta )^T \f$.
     * @param i_degreesOfFreedom of the current time step \f$ t^\text{cell} \f$ for which the time derivatives \f$ \frac{\partial^j}{\partial t^j} \f$ will be computed.
     * @param i_starMatrices star matrices, 0: \f$ A^*_k \f$, 1: \f$ B^*_k \f$, 2: \f$ C^*_k \f$.
     * @param i_timeIntegrated time integrated DOFs.
     * @param o_timeDerivatives (optional) time derivatives of the degrees of freedom in compressed format. If NULL only time integrated DOFs are returned.
     **/
    void computeAder(       double i_timeStepWidth,
                            real** i_stiffnessMatrices,
                      const real*  i_degreesOfFreedom,
                            real   i_starMatrices[3][STAR_NNZ],
                            real*  o_timeIntegrated,
                            real*  o_timeDerivatives = NULL );

    /**
     * Derives the number of non-zero and hardware floating point operation in the ADER procedure.
     * @param o_nonZeroFlops number of performed non zero floating point operations.
     * @param o_hardwareFlops number of performed floating point operations in hardware.
     **/
    void flopsAder( unsigned int &o_nonZeroFlops,
                    unsigned int &o_hardwareFlops );

    /**
     * Computes the time integrated degrees of freedom from previously computed time derivatives.
     *
     * @param i_expansionPoint expansion point (in time) of the Taylor series.
     * @param i_integrationStart start of the integration interval.
     * @param i_integrationEnd end of the integration interval.
     * @param i_timeDerivatives time derivatives.
     * @param o_timeIntegratedDofs time integrated DOFs over the interval: \f$ [ t^\text{start},  t^\text{end} ] \f$ 
     **/
    void computeIntegral( double                            i_expansionPoint,
                          double                            i_integrationStart,
                          double                            i_integrationEnd,
                          GlobalData const*,
                          seissol::model::TimeIntegrationData const*,
                          const real*                       i_timeDerivatives,
                          real                              o_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] );
                          
    
    /**
     * Convert compressed and memory aligned time derivatives to a full (including zeros) unaligned format.
     *
     * @param i_compressedDerivatives derivatives in compressed, aligned format.
     * @param o_fullDerivatives derivatives in full, unaligned format.
     **/
    template<typename real_from, typename real_to>
    static void convertAlignedCompressedTimeDerivatives( const real_from *i_compressedDerivatives,
                                                               real_to    o_fullDerivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] )
    {
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
};

#endif

