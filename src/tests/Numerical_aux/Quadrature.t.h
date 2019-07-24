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

#include <Numerical_aux/Quadrature.h>

#define EPSILON 1e-15

namespace seissol {
  namespace unit_test {
    class QuadratureTestSuite;
  }
}

class seissol::unit_test::QuadratureTestSuite : public CxxTest::TestSuite
{
public:  
	void testGaussJacobi()
	{
    double points[5];
    double weights[5];
    seissol::quadrature::GaussJacobi(points, weights, 5, 1, 3);
    // Compare to Maple reference solution
    TS_ASSERT_DELTA(points[0], 0.86698568210542769702, EPSILON);
    TS_ASSERT_DELTA(points[1], 0.57652877512667440772, EPSILON);
    TS_ASSERT_DELTA(points[2], 0.17976783188823737401, EPSILON);
    TS_ASSERT_DELTA(points[3], -.25499675973326581341, EPSILON);
    TS_ASSERT_DELTA(points[4], -.65399981510135937963, EPSILON);
    TS_ASSERT_DELTA(weights[0], 0.18915446768616357329, EPSILON);
    TS_ASSERT_DELTA(weights[1], 0.58714974961811369751, EPSILON);
    TS_ASSERT_DELTA(weights[2], 0.57657004957734461768, EPSILON);
    TS_ASSERT_DELTA(weights[3], 0.22255926867518051648, EPSILON);
    TS_ASSERT_DELTA(weights[4], 0.024566464443197594119, EPSILON);
	}

  void testGaussLobatto()
  {
    double points[5];
    double weights[5];
    seissol::quadrature::GaussLobatto(points, weights, 5);
    TS_ASSERT_DELTA(points[0], -1.0000000000000000000000000000000, EPSILON);
    TS_ASSERT_DELTA(points[1], -.65465367070797714379829245624680, EPSILON);
    TS_ASSERT_DELTA(points[2], 0.0, EPSILON);
    TS_ASSERT_DELTA(points[3], .65465367070797714379829245624685, EPSILON);
    TS_ASSERT_DELTA(points[4], 1.0000000000000000000000000000000, EPSILON);
    TS_ASSERT_DELTA(weights[0], 1.0/10.0, EPSILON);
    TS_ASSERT_DELTA(weights[1], 49.0/90.0, EPSILON);
    TS_ASSERT_DELTA(weights[2], 32.0/45.0, EPSILON);
    TS_ASSERT_DELTA(weights[3], 49.0/90.0, EPSILON);
    TS_ASSERT_DELTA(weights[4], 1.0/10.0, EPSILON);
  }
  
	void testTriangleQuadrature()
	{
    double points[4][2];
    double weights[4];
    seissol::quadrature::TriangleQuadrature(points, weights, 2);
    // Compare to Maple reference solution
    TS_ASSERT_DELTA(points[0][0], 0.64494897427831780982, EPSILON);
    TS_ASSERT_DELTA(points[1][0], 0.64494897427831780982, EPSILON);
    TS_ASSERT_DELTA(points[2][0], 0.15505102572168219018, EPSILON);
    TS_ASSERT_DELTA(points[3][0], 0.15505102572168219018, EPSILON);
    TS_ASSERT_DELTA(points[0][1], 0.28001991549907407200, EPSILON);
    TS_ASSERT_DELTA(points[1][1], 0.075031110222608118175, EPSILON);
    TS_ASSERT_DELTA(points[2][1], 0.66639024601470138669, EPSILON);
    TS_ASSERT_DELTA(points[3][1], 0.17855872826361642311, EPSILON);
    TS_ASSERT_DELTA(weights[0], 0.090979309128011415315, EPSILON);
    TS_ASSERT_DELTA(weights[1], 0.090979309128011415315, EPSILON);
    TS_ASSERT_DELTA(weights[2], 0.15902069087198858472, EPSILON);
    TS_ASSERT_DELTA(weights[3], 0.15902069087198858472, EPSILON);
	}
};
