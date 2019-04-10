/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
#include <cxxtest/TestSuite.h>

#include <Numerical_aux/Functions.h>

#define EPSILON 1e-15

namespace seissol {
  namespace unit_test {
    class FunctionsTestSuite;
  }
}

class seissol::unit_test::FunctionsTestSuite : public CxxTest::TestSuite
{
public:
	void testJacobiP()
	{
    // Compare to Maple reference solution
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 0,  1, 0,  0.5),                1.0, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 1,  0, 0, -0.3), -.3000000000000000, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 1, 13, 4, -1.0), -5.000000000000000, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP(10,  0, 4, -1.0),              1001., EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP(10,  0, 4, -0.3), -1.870307686685156, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP(10,  0, 4,  0.0), 0.8671875000000000, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP(10,  0, 4,  0.3), -.3404032083257812, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP(10,  0, 4,  1.0),                 1., EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 2,  1, 4, -1.0),                15., EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 3,  2, 3, -0.3),  1.249375000000000, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 4,  3, 2,  0.0), 0.9375000000000002, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 5,  4, 1,  0.3), -1.982012812500000, EPSILON);
    TS_ASSERT_DELTA(seissol::functions::JacobiP( 6,  5, 0,  1.0),  462.0000000000000, EPSILON);
	}
  
	void testJacobiPFirstDerivative()
	{
    double eps = 10. * EPSILON;
    // Compare to Maple reference solution
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 0,  1, 0,  0.5),                0.0, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 1,  0, 0, -0.3),  1.000000000000000, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 1, 13, 4, -0.9),  9.500000000000000, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 7,  0, 4, -0.3), -9.122014125000001, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 7,  0, 4,  0.0), 7.8750000000000, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 7,  0, 4,  0.3), -6.067144125000002, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 7,  0, 4,  0.9), -2.71497712500000, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 2,  1, 4, -0.9), -22.20000000000000, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 3,  2, 3, -0.3), 3.318750000000000, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 4,  3, 2,  0.0), -3.749999999999998, eps);
    TS_ASSERT_DELTA(seissol::functions::JacobiPFirstDerivative( 5,  4, 1,  0.3), -15.37232812500000, eps);
	}

  void testTetraDubinerP()
  {
    std::vector<std::vector<double>> tests = {
      {0, 0, 0, 0.25, 0.25, 0.25, 1},
      {1, 0, 0, 0.25, 0.25, 0.25, 0.},
      {0, 1, 0, 0.25, 0.25, 0.25, 0.},
      {0, 0, 1, 0.25, 0.25, 0.25, 0.},
      {2, 0, 0, 0.25, 0.25, 0.25, -0.1250},
      {1, 1, 0, 0.25, 0.25, 0.25, 0.},
      {0, 2, 0, 0.25, 0.25, 0.25, -0.3125},
      {1, 0, 1, 0.25, 0.25, 0.25, 0.},
      {0, 1, 1, 0.25, 0.25, 0.25, 0.},
      {0, 0, 2, 0.25, 0.25, 0.25, -0.5625},
      {3, 0, 0, 0.25, 0.25, 0.25, -0.},
      {2, 1, 0, 0.25, 0.25, 0.25, -0.125000},
      {1, 2, 0, 0.25, 0.25, 0.25, -0.},
      {0, 3, 0, 0.25, 0.25, 0.25, 0.125000},
      {2, 0, 1, 0.25, 0.25, 0.25, -0.125000},
      {1, 1, 1, 0.25, 0.25, 0.25, 0.},
      {0, 2, 1, 0.25, 0.25, 0.25, -0.31250000000000000000},
      {1, 0, 2, 0.25, 0.25, 0.25, -0.},
      {0, 1, 2, 0.25, 0.25, 0.25, -0.},
      {0, 0, 3, 0.25, 0.25, 0.25, 0.437500},
      {0, 0, 0, 0., 0.5, 0.5, 1.},
      {1, 0, 0, 0., 0.5, 0.5, 0.},
      {0, 1, 0, 0., 0.5, 0.5, 1.0},
      {0, 0, 1, 0., 0.5, 0.5, 1.0},
      {2, 0, 0, 0., 0.5, 0.5, 0.},
      {1, 1, 0, 0., 0.5, 0.5, 0.},
      {0, 2, 0, 0., 0.5, 0.5, 0.75},
      {1, 0, 1, 0., 0.5, 0.5, 0.},
      {0, 1, 1, 0., 0.5, 0.5, 2.00},
      {0, 0, 2, 0., 0.5, 0.5, -0.25},
      {3, 0, 0, 0., 0.5, 0.5, 0.},
      {2, 1, 0, 0., 0.5, 0.5, 0.},
      {1, 2, 0, 0., 0.5, 0.5, 0.},
      {0, 3, 0, 0., 0.5, 0.5, 0.500},
      {2, 0, 1, 0., 0.5, 0.5, 0.},
      {1, 1, 1, 0., 0.5, 0.5, 0.},
      {0, 2, 1, 0., 0.5, 0.5, 2.2500000000000000000},
      {1, 0, 2, 0., 0.5, 0.5, 0.},
      {0, 1, 2, 0., 0.5, 0.5, 1.0},
      {0, 0, 3, 0., 0.5, 0.5, -0.750},
      {0, 0, 0, 0., 1.0, 0., 1.},
      {1, 0, 0, 0., 1.0, 0., 0.},
      {0, 1, 0, 0., 1.0, 0., 2.0},
      {0, 0, 1, 0., 1.0, 0., -1.},
      {2, 0, 0, 0., 1.0, 0., 0.},
      {1, 1, 0, 0., 1.0, 0., 0.},
      {0, 2, 0, 0., 1.0, 0., 3.00},
      {1, 0, 1, 0., 1.0, 0., -0.},
      {0, 1, 1, 0., 1.0, 0., -2.0},
      {0, 0, 2, 0., 1.0, 0., 1.},
      {3, 0, 0, 0., 1.0, 0., 0.},
      {2, 1, 0, 0., 1.0, 0., 0.},
      {1, 2, 0, 0., 1.0, 0., 0.},
      {0, 3, 0, 0., 1.0, 0., 4.000},
      {2, 0, 1, 0., 1.0, 0., -0.},
      {1, 1, 1, 0., 1.0, 0., -0.},
      {0, 2, 1, 0., 1.0, 0., -3.0000000000000000000},
      {1, 0, 2, 0., 1.0, 0., 0.},
      {0, 1, 2, 0., 1.0, 0., 2.0},
      {0, 0, 3, 0., 1.0, 0., -1.}
    };
    for (auto const& t : tests) {
      TS_TRACE(t);
      TS_ASSERT_DELTA(seissol::functions::TetraDubinerP(static_cast<unsigned>(t[0]), static_cast<unsigned>(t[1]), static_cast<unsigned>(t[2]), t[3], t[4], t[5]), t[6], EPSILON);
    }
  }
};
