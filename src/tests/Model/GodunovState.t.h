/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 * @section LICENSE
 * Copyright (c) 2015-2019, SeisSol Group
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
#include <cassert>
#include <cmath>

#include <Model/common.hpp>
#include <Model/Setup.h>
#include <Model/common_datastructures.hpp>


real const EPSILON = std::numeric_limits<real>::epsilon();

namespace seissol {
  namespace unit_test {
    class GodunovStateTestSuite;
  }
}

class seissol::unit_test::GodunovStateTestSuite : public CxxTest::TestSuite
{
  public:
    void testGodunovState()
    {
      seissol::model::Material local;
      seissol::model::Material neighbor;

      real localData[seissol::tensor::QgodLocal::size()];
      real neighborData[seissol::tensor::QgodLocal::size()];
      init::QgodLocal::view::type QgodLocal = init::QgodLocal::view::create(localData);
      init::QgodNeighbor::view::type QgodNeighbor = init::QgodNeighbor::view::create(neighborData); 
      QgodLocal.setZero();
      QgodNeighbor.setZero();

      //test homogeneous material
      local.rho = 2700;
      local.mu = 3.23980992e10;
      local.lambda = 3.24038016e10;
      neighbor.rho = 2700;
      neighbor.mu = 3.23980992e10;
      neighbor.lambda = 3.24038016e10;

      std::array<std::array<double, 9>, 9> solution_homogeneous = {{
        {  0.5000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  8.100000000000001e6,   0.0000000000000000,   0.0000000000000000},
        {  0.1666862222222223,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  2.700316800000001e6,   0.0000000000000000,   0.0000000000000000},
        {  0.1666862222222223,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  2.700316800000001e6,   0.0000000000000000,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.5000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  4.676399999999999e6,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.5000000000000000,   0.0000000000000000,   0.0000000000000000,  4.676399999999999e6},
        {3.086419753086419e-8,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.5000000000000000,   0.0000000000000000,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000, 5.345992643914123e-8,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.4999999999999999,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000, 5.345992643914123e-8,   0.0000000000000000,   0.0000000000000000,   0.4999999999999999}
      }};

      seissol::model::getTransposedGodunovState( local,
          neighbor,
          regular,
          QgodLocal,
          QgodNeighbor );
      test_local(QgodLocal, solution_homogeneous);
      test_neighbor(QgodNeighbor, solution_homogeneous);

      //test free surface
      std::array<std::array<double, 9>, 9> solution_boundary = {{
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000},
        {-6.172839506172840e-8,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,     1.000000000000000,    0.0000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000, -1.069198528782824e-7,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000,     1.000000000000000,    0.0000000000000000},
        {   0.0000000000000000,   0.0000000000000000,    0.0000000000000000,    0.0000000000000000,    0.0000000000000000, -1.069198528782824e-7,    0.0000000000000000,    0.0000000000000000,     1.000000000000000}
      }};
      seissol::model::getTransposedGodunovState( local,
          neighbor,
          freeSurface,
          QgodLocal,
          QgodNeighbor );
      test_neighbor(QgodLocal, solution_boundary);
      test_nan(QgodNeighbor);

      
      
      //test heterogeneous material
      local.rho = 2700;
      local.mu = 3.23980992e10;
      local.lambda = 3.24038016e10;
      neighbor.rho = 2600;
      neighbor.mu = 1.04e10;
      neighbor.lambda = 2.08e10;

      std::array<std::array<double, 9>, 9> solution_heterogeneous = {{
        {  0.6090225563909775,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  6.333834586466166e6,   0.0000000000000000,   0.0000000000000000},
        {  0.2030313383458647,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  2.111525918796993e6,   0.0000000000000000,   0.0000000000000000},
        {  0.2030313383458647,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  2.111525918796993e6,   0.0000000000000000,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.6426804463745808,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,  3.341938321147820e6,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.6426804463745808,   0.0000000000000000,   0.0000000000000000,  3.341938321147820e6},
        {3.759398496240601e-8,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.3909774436090225,   0.0000000000000000,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000, 6.871529877411907e-8,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.3573195536254192,   0.0000000000000000},
        {  0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000,   0.0000000000000000, 6.871529877411907e-8,   0.0000000000000000,   0.0000000000000000,   0.3573195536254192}
      }};

      seissol::model::getTransposedGodunovState( local,
          neighbor,
          regular,
          QgodLocal,
          QgodNeighbor );
      test_local(QgodLocal, solution_heterogeneous);
      test_neighbor(QgodNeighbor, solution_heterogeneous);

    }



    void test_local(init::QgodLocal::view::type QgodLocal, std::array<std::array<double, 9>, 9> solution) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          if (i == j) {
            auto abs_solution = std::abs(1-solution[j][i]);
            auto abs_calculated = std::abs(QgodLocal(i,j));
            TS_ASSERT_DELTA(QgodLocal(i,j), 1-solution[j][i], abs_solution*1e2*EPSILON);
            TS_ASSERT_DELTA(QgodLocal(i,j), 1-solution[j][i], abs_calculated*1e2*EPSILON);
          }
          else {
            auto abs_solution = std::abs(solution[j][i]);
            auto abs_calculated = std::abs(QgodLocal(i,j));
            TS_ASSERT_DELTA(QgodLocal(i,j), -solution[j][i], abs_solution*1e2*EPSILON);
            TS_ASSERT_DELTA(QgodLocal(i,j), -solution[j][i], abs_calculated*1e2*EPSILON);
          }
        }
      } 
    }

    void test_neighbor(init::QgodLocal::view::type QgodNeighbor, std::array<std::array<double, 9>, 9> solution) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            auto abs_solution = std::abs(solution[j][i]);
            auto abs_calculated = std::abs(QgodNeighbor(i,j));
            TS_ASSERT_DELTA(QgodNeighbor(i,j), solution[j][i], abs_solution*1e2*EPSILON);
            TS_ASSERT_DELTA(QgodNeighbor(i,j), solution[j][i], abs_calculated*1e2*EPSILON);
        }
      }
    }

    void test_nan(init::QgodLocal::view::type QgodNeighbor) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          //CXXTEST TS_ASSERT_IS_NAN is broken
          assert(std::isnan(QgodNeighbor(i,j)));
        }
      }
    }
};
