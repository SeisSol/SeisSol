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
     * Gets the lts setup in relation to the four face neighbors.
     *   Remark: Remember to perform the required normalization step.
     *
     * -------------------------------------------------------------------------------
     *
     *  0 in one of the first four bits: Face neighboring data are buffers.
     *  1 in one of the first four bits: Face neighboring data are derivatives.
     *
     *     Example 1:
     *     [           12 rem. bits               | buf/der bits ]
     *     [  -  -  -  -  -  -  -  -  -  -  -  -  |  0  1  1  0  ]
     *     [ 15 14 13 12 11 10  9  8  7  6  5  4  |  3  2  1  0  ]
     *  In Example 1 the data for face neighbors 0 and 3 are buffers and for 1 and 2 derivatives.
     *
     *  0 in one of bits 4 - 7: No global time stepping
     *  1 in one of bits 4 - 7: The  current cell has a global time stepping relation with the face neighbor.
     *
     *     Example 2:
     *     [       8 rem. bits       |   GTS bits  | buf/der bits ]
     *     [  -  -  -  -  -  -  -  - | 0  0  1  1  |  0  1  1  0  ]
     *     [ 15 14 13 12 11 10  9  8 | 7  6  5  4  |  3  2  1  0  ]
     *  In Example 2 the data of face neighbors 0 and 3 are buffers, 1 and 2 deliver derivatives
     *  Face neighbor 0 has a GTS-relation and this cell works directly on the delivered buffer.
     *  Face neighbor 1 has a GTS-relation, but delivers derivatives -> The derivatives have to translated to time integrated DOFs first.
     *  Face neighbor 2 has a LTS-relation and receives derivatives from its neighbor -> The derivates have to be used for a partial time integration.
     *  Face neighbor 3 has a LTS-relation and can operate on the buffers directly.
     *
     * -------------------------------------------------------------------------------
     *
     *  1 in the eigth bit: the cell is required to work on time integration buffers.
     *  1 in the nineth bit: the cell is required to compute time derivatives.
     *
     *     Example 3:
     *     [     remaining     | der. buf. |       first 8 bits       ]
     *     [  -  -  -  -  -  - |  0    1   |  -  -  -  -  -  -  -  -  ]
     *     [ 15 14 13 12 11 10 |  9    8   |  7  6  5  4  3  2  1  0  ]
     *  In Example 3 only a buffer is stored as for example in global time stepping.
     *
     *     Example 4:
     *     [     remaining     | der. buf. |       first 8 bits       ]
     *     [  -  -  -  -  -  - |  1    1   |  -  -  -  -  -  -  -  -  ]
     *     [ 15 14 13 12 11 10 |  9    8   |  7  6  5  4  3  2  1  0  ]
     *  In Example 4 both (buffer+derivative) is stored.
     *
     * -------------------------------------------------------------------------------
     *
     *  1 in the tenth bit: the cell local buffer is a LTS buffer (reset on request only).
     *
     *     Example 5:
     *     [   remaining    | LTS buf. |          first 10 bits         ]
     *     [  -  -  -  -  - |     1    |  -  1  -  -  -  -  -  -  -  -  ]
     *     [ 15 14 13 12 11 |    10    |  9  8  7  6  5  4  3  2  1  0  ]
     *  In Example 5 the buffer is a LTS buffer (reset on request only). GTS buffers are updated in every time step.
     *
     * @return lts setup.
     * @param i_localCluster global id of the cluster to which this cell belongs.
     * @param i_neighboringClusterIds global ids of the clusters the face neighbors belong to (if present).
     * @param i_faceTypes types of the four faces.
     * @param i_faceNeighborIds face neighbor ids.
     * @param i_copy true if the cell is part of the copy layer (only required for correctness in dynamic rupture computations).
     **/
    static unsigned short getLtsSetup(      unsigned int   i_localClusterId,
                                            unsigned int   i_neighboringClusterIds[4],
                                      const enum faceType  i_faceTypes[4],
                                      const unsigned int   i_faceNeighborIds[4], // TODO: Remove, outdated
                                            bool           i_copy = false ) {
      // reset the LTS setup
      unsigned short l_ltsSetup = 0;

      // iterate over the faces
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // continue for boundary conditions
        if( i_faceTypes[l_face] != regular        &&
            i_faceTypes[l_face] != dynamicRupture &&
            i_faceTypes[l_face] != periodic       &&
            i_faceTypes[l_face] != freeSurface ) {
          continue;
        }
        // free surface fake neighbors are GTS
        else if( i_faceTypes[l_face] == freeSurface ) {
          l_ltsSetup |= (1 << (l_face+4) );
        }
        // dynamic rupture faces are always global time stepping but operate on derivatives
        else if( i_faceTypes[l_face] == dynamicRupture ) {
          // face-neighbor provides GTS derivatives
          // face-neighbor provides derivatives
          l_ltsSetup |= ( 1 <<  l_face      );
          l_ltsSetup |= ( 1 << (l_face + 4) );

          // cell is required to provide derivatives for dynamic rupture
          l_ltsSetup |= ( 1 << 9 );

          if( i_copy ) { // set the buffer invalid in copy layers
                         // TODO: Minor improvements possible: Non-DR MPI-neighbor for example
            l_ltsSetup |= ( 1 << 10 );
          }
        }
        // derive the LTS setup based on the cluster ids
        else {
          // neighboring cluster has a larger time step than this cluster
          if( i_localClusterId < i_neighboringClusterIds[l_face] ) {
            // neighbor delivers time derivatives
            l_ltsSetup |= ( 1 << l_face );

            // the cell-local buffer is used in LTS-fashion
            l_ltsSetup |= ( 1 << 10     );
          }
          // GTS relation
          else if( i_localClusterId == i_neighboringClusterIds[l_face] ) {
            l_ltsSetup |= ( 1 << (l_face + 4) );
          }

          // cell is required to provide derivatives
          if( i_localClusterId > i_neighboringClusterIds[l_face] ) {
            l_ltsSetup |= ( 1 << 9 );
          }
          // cell is required to provide a buffer
          else {
            l_ltsSetup |= ( 1 << 8 );
          }
        }

        // true lts buffer with gts required derivatives
        if( (l_ltsSetup >> 10)%2 == 1 && (l_ltsSetup >> 4)%16 != 0 ) {
          l_ltsSetup |= ( 1 << 9 );
        }
      }

      /*
       * Normalize for special case "free surface on derivatives":
       *   If a cell provides either buffers in a LTS fashion or derivatives only,
       *   the neighboring contribution of the boundary intergral is required to work on the cells derivatives.
       *   Free surface boundary conditions work on the cells DOFs in the neighboring contribution:
       *   Enable cell local derivatives in this case and mark that the "fake neighbor" provides derivatives.
       */
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // check for special case free-surface requirements
        if( i_faceTypes[l_face] == freeSurface &&       // free surface face
           ( (l_ltsSetup >> 10) % 2 == 1 ||             // lts fashion buffer
             (l_ltsSetup >> 8 ) % 2 == 0 )         ) {  // no buffer at all
          l_ltsSetup |= ( 1 << 9 );       // enable derivatives computation
          l_ltsSetup |= ( 1 << l_face );  // enable derivatives for the fake face neighbor
        }
      }

      return l_ltsSetup;
    }

    /**
     * Normalizes the LTS setup for the special case "GTS on derivatives":
     *   If a face neighbor provides true buffers to cells with larger time steps,
     *   the local cell is required to operate on derivatives of this face neighbor.
     *
     *   Example:
     *        | own |  fn 1 |  fn 2 | fn 3 | fn 4 |
     *   local|  dt |    dt | 0.5dt |   dt |   dt |
     *   fn 4 |  dt | 0.5dt |   2dt |   dt |   dt |
     *         -----------------------------------
     *   In the example the local cell is connected via face 4 with to a GTS neighbor.
     *   Face neighbor 4 is required to deliver true buffers to its second face neighbor.
     *   It follows that the local cell has to operate on the derivatives of face neighbor 4.
     *
     * @param i_neighboringSetups local time stepping setups for the neighboring cells, set to GTS (240) if not defined (e.g. in case of boundary conditions).
     * @param io_localLtsSetup local time stepping setup of the local cell.
     **/
    static void normalizeLtsSetup( unsigned short  i_neighboringLtsSetups[4],
                                   unsigned short &io_localLtsSetup ) {
      // iterate over the face neighbors
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // enforce derivatives if this is a "GTS on derivatives" relation
        if( (io_localLtsSetup >> (l_face + 4))%2 && (i_neighboringLtsSetups[l_face] >> 10)%2 == 1 ) {
          io_localLtsSetup |= (1 << l_face);
        }
      }
    }

    /**
     * Gets the total size of the time derivative for a specific convergence order and alignment.
     *
     * @param i_convergenceOrder convergence order of the ADER-DG method.
     * @param i_numberOfQuantities number of quantities.
     * @param i_alignment memory alignment in bits.
     **/
    static unsigned int getAlignedTimeDerivativesSize( unsigned int i_convergenceOrder   = CONVERGENCE_ORDER,
                                                       unsigned int i_numberOfQuantities = NUMBER_OF_QUANTITIES,
                                                       unsigned int i_alignment          = ALIGNMENT ) {
      assert( i_convergenceOrder >= 1 );
      unsigned int l_alignedDerivativesSize = 0;

      for( unsigned int l_associatedOrder = i_convergenceOrder; l_associatedOrder > 0; l_associatedOrder-- ) {
        l_alignedDerivativesSize += getNumberOfAlignedBasisFunctions( l_associatedOrder, i_alignment ) * i_numberOfQuantities;
      }

      return l_alignedDerivativesSize;
    }

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
     * Evaluates the taylor series expansion at the given evaluation point.
     *
     * @param i_expansionPoint expansion point of the taylor series (point in time where the time derivatives were computed from).
     * @param i_evaluationPoint point in time when the taylor series is evaluated.
     * @param i_timeDerivatives time derivatives.
     * @param o_degreesOfFreedom degrees of freedom at the time of the evaluation point.
     **/
    void computeExtrapolation(       double i_expansionPoint,
                                     double i_evaluationPoint,
                               const real*  i_timeDerivatives,
                                     real*  o_degreesOfFreedom );

    /**
     * Computes the time integrated degrees of freedom from previously computed time derivatives.
     *
     * @param i_expansionPoint expansion point (in time) of the Taylor series.
     * @param i_integrationStart start of the integration interval.
     * @param i_integrationEnd end of the integration interval.
     * @param i_timeDerivatives time derivatives.
     * @param o_timeIntegratedDofs time integrated DOFs over the interval: \f$ [ t^\text{start},  t^\text{end} ] \f$ 
     **/
    void computeIntegral(       double i_expansionPoint,
                                double i_integrationStart,
                                double i_integrationEnd,
                          const real*  i_timeDerivatives,
                                real   o_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] );

    /**
     * Either copies pointers to the DOFs in the time buffer or integrates the DOFs via time derivatives.
     *   Evaluation depends on bit 0-3  of the LTS setup.
     *   0 -> copy buffer; 1 -> integrate via time derivatives
     *     Example:

     *     [     4 unused     | copy or int bits  ]
     *     [ -    -    -    - |  0    1    1    0 ]
     *     [ 7    6    5    4 |  3    2    1    0 ]
     *
     *   0 - 0: time integrated DOFs of cell 0 are copied from the buffer.
     *   1 - 1: DOFs of cell 1 are integrated in time via time derivatives.
     *   2 - 1: DOFs of cell 2 are integrated in time via time derivaitves.
     *   3 - 0: time itnegrated DOFs of cell 3 are copied from the buffer.
     *
     * @param i_ltsSetup bitmask for the LTS setup.
     * @param i_faceTypes face types of the neighboring cells.
     * @param i_currentTime current time of the cell [0] and it's four neighbors [1], [2], [3] and [4].
     * @param i_timeStepWidth time step width of the cell.
     * @param i_timeDofs pointers to time integrated buffers or time derivatives of the four neighboring cells.
     * @param i_integrationBuffer memory where the time integration goes if derived from derivatives. Ensure thread safety!
     * @param o_timeIntegrated pointers to the time integrated DOFs of the four neighboring cells (either local integration buffer or integration buffer of input).
     **/
    void computeIntegrals( unsigned short      i_ltsSetup,
                           const enum faceType i_faceTypes[4],
                           const double        i_currentTime[5],
                           double              i_timeStepWidth,
                           real * const        i_timeDofs[4],
                           real                o_integrationBuffer[4][NUMBER_OF_ALIGNED_DOFS],
                           real *              o_timeIntegrated[4] );

    /**
     * Special case of the computeIntergals function, which assumes a common "current time" for all face neighbors which provide derivatives.
     *
     * @param i_ltsSetup bitmask for the LTS setup.
     * @param i_faceTypes face types of the neighboring cells.
     * @param i_timeStepStart start time of the current cell with respect to the common point zero: Time of the larger time step width prediction of the face neighbors.
     * @param i_timeStepWidth time step width of the cell.
     * @param i_timeDofs pointers to time integrated buffers or time derivatives of the four neighboring cells.
     * @param i_integrationBuffer memory where the time integration goes if derived from derivatives. Ensure thread safety!
     * @param o_timeIntegrated pointers to the time integrated DOFs of the four neighboring cells (either local integration buffer or integration buffer of input).
     **/
    void computeIntegrals( unsigned short      i_ltsSetup,
                           const enum faceType i_faceTypes[4],
                           const double        i_timeStepStart,
                           const double        i_timeStepWidth,
                           real * const        i_timeDofs[4],
                           real                o_integrationBuffer[4][NUMBER_OF_ALIGNED_DOFS],
                           real *              o_timeIntegrated[4] );
};

#endif

