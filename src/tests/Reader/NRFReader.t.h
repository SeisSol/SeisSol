/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 **/

#include <iostream>

#include <cxxtest/TestSuite.h>

#include <SourceTerm/NRFReader.h>
#include <SourceTerm/NRF.h>

#include "slipRates.h"

#if defined(DOUBLE_PRECISION)
#define EPSILON 1e-16
#elif defined(SINGLE_PRECISION)
#define EPSILON 1e-8
#endif

namespace seissol {
  namespace unit_test {
    class NRFReaderTestSuite;
  }
}

class seissol::unit_test::NRFReaderTestSuite : public CxxTest::TestSuite
{
public:
	void testTransformMomentTensor()
	{
          seissol::sourceterm::NRF nrf;
          seissol::sourceterm::readNRF("Testing/source_loh.nrf", nrf);

          TS_ASSERT_DELTA(nrf.centres[0].x,    0.0, EPSILON);
          TS_ASSERT_DELTA(nrf.centres[0].y,    0.0, EPSILON);
          TS_ASSERT_DELTA(nrf.centres[0].z, 2000.0, EPSILON);
          
          TS_ASSERT_DELTA(nrf.subfaults[0].tan1.x, 0.0, EPSILON); 
          TS_ASSERT_DELTA(nrf.subfaults[0].tan1.y, 1.0, EPSILON); 
          TS_ASSERT_DELTA(nrf.subfaults[0].tan1.z, 0.0, EPSILON); 

          TS_ASSERT_DELTA(nrf.subfaults[0].tan2.x, 0.0, EPSILON); 
          TS_ASSERT_DELTA(nrf.subfaults[0].tan2.y, 0.0, EPSILON); 
          TS_ASSERT_DELTA(nrf.subfaults[0].tan2.z, 1.0, EPSILON); 

          TS_ASSERT_DELTA(nrf.subfaults[0].normal.x, 1.0, EPSILON); 
          TS_ASSERT_DELTA(nrf.subfaults[0].normal.y, 0.0, EPSILON); 
          TS_ASSERT_DELTA(nrf.subfaults[0].normal.z, 0.0, EPSILON);

          TS_ASSERT_DELTA(nrf.subfaults[0].area, 3.0866008336686616479e+07, EPSILON);
          TS_ASSERT_DELTA(nrf.subfaults[0].tinit, 0.0, EPSILON);
          TS_ASSERT_DELTA(nrf.subfaults[0].timestep, 0.0002, EPSILON);
          TS_ASSERT_DELTA(nrf.subfaults[0].mu, 0.0, EPSILON);
          TS_ASSERT_EQUALS(nrf.source, 1);

	  for (int i = 0; i < 80; i++) {
	    //NRF File is in cm
	    TS_ASSERT_DELTA(nrf.sliprates[0][i], slipRates[i]/100, EPSILON);
	    TS_ASSERT_DELTA(nrf.sliprates[1][i], 0, EPSILON);
	    TS_ASSERT_DELTA(nrf.sliprates[2][i], 0, EPSILON);
	  }


        }
};
