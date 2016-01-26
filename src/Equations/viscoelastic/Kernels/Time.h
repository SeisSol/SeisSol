/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
#include <generated_code/sizes.h>

namespace seissol {
  namespace kernels {
    class Time;
  }
}

class seissol::kernels::Time {
  public:
    /**
     * Constructor, which initializes the time kernel.
     **/
    Time() {}

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
    void computeAder( double                i_timeStepWidth,
                      GlobalData*           global,
                      LocalIntegrationData* local,
                      real const*           i_degreesOfFreedom,
                      real*                 o_timeIntegrated,
                      real*                 o_timeDerivatives = NULL );

    /**
     * Derives the number of non-zero and hardware floating point operation in the ADER procedure.
     * @param o_nonZeroFlops number of performed non zero floating point operations.
     * @param o_hardwareFlops number of performed floating point operations in hardware.
     **/
    void flopsAder( unsigned int &o_nonZeroFlops,
                    unsigned int &o_hardwareFlops );
                    
    unsigned bytesAder();

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
     * Convert compressed and memory aligned time derivatives to a full (including zeros) unaligned format.
     *
     * @param i_compressedDerivatives derivatives in compressed, aligned format.
     * @param o_fullDerivatives derivatives in full, unaligned format.
     **/
    template<typename real_from, typename real_to>
    static void convertAlignedCompressedTimeDerivatives( const real_from *i_compressedDerivatives,
                                                               real_to    o_fullDerivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] )
    {
        for (unsigned order = 0; order < CONVERGENCE_ORDER; ++order) {
          seissol::kernels::copySubMatrix( &i_compressedDerivatives[order * NUMBER_OF_ALIGNED_DOFS],
                                           NUMBER_OF_BASIS_FUNCTIONS,
                                           NUMBER_OF_QUANTITIES,
                                           NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                           o_fullDerivatives[order],
                                           NUMBER_OF_BASIS_FUNCTIONS,
                                           NUMBER_OF_QUANTITIES,
                                           NUMBER_OF_BASIS_FUNCTIONS );
        }
    }
};

#endif

