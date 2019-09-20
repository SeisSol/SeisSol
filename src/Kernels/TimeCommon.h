/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 * Common functions for SeisSol's time kernel.
 **/

#ifndef KERNELS_TIMECOMMON_H_
#define KERNELS_TIMECOMMON_H_

#include <Initializer/typedefs.hpp>
#include <Kernels/Time.h>
#include <generated_code/tensor.h>

namespace seissol {
  namespace kernels {
    namespace TimeCommon {
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
      void computeIntegrals(  Time&                             i_time,
                              unsigned short                    i_ltsSetup,
                              const enum faceType               i_faceTypes[4],
                              const double                      i_currentTime[5],
                              double                            i_timeStepWidth,
                              real * const                      i_timeDofs[4],
                              real                              o_integrationBuffer[4][tensor::I::size()],
                              real *                            o_timeIntegrated[4] );

      /**
       * Special case of the computeIntegrals function, which assumes a common "current time" for all face neighbors which provide derivatives.
       *
       * @param i_ltsSetup bitmask for the LTS setup.
       * @param i_faceTypes face types of the neighboring cells.
       * @param i_timeStepStart start time of the current cell with respect to the common point zero: Time of the larger time step width prediction of the face neighbors.
       * @param i_timeStepWidth time step width of the cell.
       * @param i_timeDofs pointers to time integrated buffers or time derivatives of the four neighboring cells.
       * @param i_integrationBuffer memory where the time integration goes if derived from derivatives. Ensure thread safety!
       * @param o_timeIntegrated pointers to the time integrated DOFs of the four neighboring cells (either local integration buffer or integration buffer of input).
       **/
      void computeIntegrals(  Time&                             i_time,
                              unsigned short                    i_ltsSetup,
                              const enum faceType               i_faceTypes[4],
                              const double                      i_timeStepStart,
                              const double                      i_timeStepWidth,
                              real * const                      i_timeDofs[4],
                              real                              o_integrationBuffer[4][tensor::I::size()],
                              real *                            o_timeIntegrated[4] );
    }
  }
}

#endif

