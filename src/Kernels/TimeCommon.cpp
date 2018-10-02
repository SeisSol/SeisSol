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
 **/

#include "TimeCommon.h"
#include <stdint.h>

void seissol::kernels::TimeCommon::computeIntegrals(  Time&                             i_time,
                                                      unsigned short                    i_ltsSetup,
                                                      const enum faceType               i_faceTypes[4],
                                                      const double                      i_currentTime[5],
                                                      double                            i_timeStepWidth,
                                                      real * const                      i_timeDofs[4],
                                                      real                              o_integrationBuffer[4][tensor::I::size()],
                                                      real *                            o_timeIntegrated[4] )
{
  /*
   * assert valid input.
   */
  // only lower 10 bits are used for lts encoding
  assert (i_ltsSetup < 2048 );

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for( int l_dofeighbor = 0; l_dofeighbor < 4; l_dofeighbor++ ) {
    assert( ((uintptr_t)i_timeDofs[l_dofeighbor])          % ALIGNMENT == 0 );
    assert( ((uintptr_t)o_integrationBuffer[l_dofeighbor]) % ALIGNMENT == 0 );
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for( unsigned int l_dofeighbor = 0; l_dofeighbor < 4; l_dofeighbor++ ) {
    // collect information only in the case that neighboring element contributions are required
    if( i_faceTypes[l_dofeighbor] != outflow && i_faceTypes[l_dofeighbor] != dynamicRupture ) {
      // check if the time integration is already done (-> copy pointer)
      if( (i_ltsSetup >> l_dofeighbor ) % 2 == 0 ) {
        o_timeIntegrated[l_dofeighbor] = i_timeDofs[l_dofeighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        i_time.computeIntegral( i_currentTime[    l_dofeighbor+1],
                                i_currentTime[    0           ],
                                i_currentTime[    0           ] + i_timeStepWidth,
                                i_timeDofs[       l_dofeighbor],
                                o_integrationBuffer[ l_dofeighbor] );

        o_timeIntegrated[l_dofeighbor] = o_integrationBuffer[ l_dofeighbor];
      }
    }
  }
}

void seissol::kernels::TimeCommon::computeIntegrals(  Time&                             i_time,
                                                      unsigned short                    i_ltsSetup,
                                                      const enum faceType               i_faceTypes[4],
                                                      const double                      i_timeStepStart,
                                                      const double                      i_timeStepWidth,
                                                      real * const                      i_timeDofs[4],
                                                      real                              o_integrationBuffer[4][tensor::I::size()],
                                                      real *                            o_timeIntegrated[4] )
{
  double l_startTimes[5];
  l_startTimes[0] = i_timeStepStart;
  l_startTimes[1] = l_startTimes[2] = l_startTimes[3] = l_startTimes[4] = 0;

  // adjust start times for GTS on derivatives
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    if( (i_ltsSetup >> (l_face + 4) ) % 2 ) {
      l_startTimes[l_face+1] = i_timeStepStart;
    }
  }

  // call the more general assembly
  computeIntegrals( i_time,
                    i_ltsSetup,
                    i_faceTypes,
                    l_startTimes,
                    i_timeStepWidth,
                    i_timeDofs,
                    o_integrationBuffer,
                    o_timeIntegrated );
}
