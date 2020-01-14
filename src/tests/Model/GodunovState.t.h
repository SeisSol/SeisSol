#include <cxxtest/TestSuite.h>
#include <cassert>
#include <cmath>

#include <Model/common_datastructures.hpp>
#include <Model/common.hpp>
#include <Model/Setup.h>

namespace seissol {
  namespace unit_test {
    class GodunovStateTestSuite;
  }
}

class seissol::unit_test::GodunovStateTestSuite : public CxxTest::TestSuite
{
  public:
    const real epsilon = std::numeric_limits<real>::epsilon();

    void testGodunovState()
    {
      real localData[seissol::tensor::QgodLocal::size()];
      real neighborData[seissol::tensor::QgodLocal::size()];
      init::QgodLocal::view::type QgodLocal = init::QgodLocal::view::create(localData);
      init::QgodNeighbor::view::type QgodNeighbor = init::QgodNeighbor::view::create(neighborData); 
      QgodLocal.setZero();
      QgodNeighbor.setZero();

      //test homogeneous material
#ifdef USE_ANISOTROPIC
      seissol::model::AnisotropicMaterial local;
      seissol::model::AnisotropicMaterial neighbor;
      double materialVal_1[22] = {
        2700,
        97200000000,
        32403801600,
        32403801600,
        0,
        0,
        0,
        97200000000,
        32403801600,
        0,
        0,
        0,
        97200000000,
        0,
        0,
        0,
        32398099200,
        0,
        0,
        32398099200,
        0,
        32398099200
      };
      setMaterial(materialVal_1, 22, &local);
      setMaterial(materialVal_1, 22, &neighbor);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
      seissol::model::ViscoElasticMaterial local;
      seissol::model::ViscoElasticMaterial neighbor;
      double materialVal_1[3 + NUMBER_OF_RELAXATION_MECHANISMS * 4];
      materialVal_1[0] = 2700;
      materialVal_1[1] = 3.23980992e10;
      materialVal_1[2] = 3.24038016e10;
      for (int i = 3; i < 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4; i++) {
        materialVal_1[i] = 0;  
      }
      setMaterial(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4, &local);
      setMaterial(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4, &neighbor);
#else
      seissol::model::ElasticMaterial local;
      seissol::model::ElasticMaterial neighbor;
      double materialVal_1[3] = {
        2700,
        3.23980992e10,
        3.24038016e10
      };
      setMaterial(materialVal_1, 3, &local);
      setMaterial(materialVal_1, 3, &neighbor);
#endif

      std::array<std::array<real, 9>, 9> solution_homogeneous = {{
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
          FaceType::regular,
          QgodLocal,
          QgodNeighbor );
      test_local(QgodLocal, solution_homogeneous);
      test_neighbor(QgodNeighbor, solution_homogeneous);

      //test free surface
      std::array<std::array<real, 9>, 9> solution_boundary = {{
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
          FaceType::freeSurface,
          QgodLocal,
          QgodNeighbor );
      test_neighbor(QgodLocal, solution_boundary);
      test_nan(QgodNeighbor);



      //test heterogeneous material
#ifdef USE_ANISOTROPIC
      double materialVal_2[22] = { 
        2600.0,
        41600000000,
        20800000000,
        20800000000,
        0,
        0,
        0,
        41600000000,
        20800000000,
        0,
        0,
        0,
        41600000000,
        0,
        0,
        0,
        10400000000,
        0,
        0,
        10400000000,
        0,
        10400000000
      };
      setMaterial(materialVal_2, 22, &neighbor);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
      double materialVal_2[3 + NUMBER_OF_RELAXATION_MECHANISMS * 4];
      materialVal_2[0] = 2600;
      materialVal_2[1] = 1.04e10;
      materialVal_2[2] = 2.08e10;
      for (int i = 3; i < 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4; i++) {
        materialVal_2[i] = 0;  
      }
      setMaterial(materialVal_2, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4, &neighbor);
#else
      double materialVal_2[3] = {
        2600,
        1.04e10,
        2.08e10
      };
      setMaterial(materialVal_2, 3, &neighbor);
#endif

      std::array<std::array<real, 9>, 9> solution_heterogeneous = {{
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
          FaceType::regular,
          QgodLocal,
          QgodNeighbor );
      test_local(QgodLocal, solution_heterogeneous);
      test_neighbor(QgodNeighbor, solution_heterogeneous);

    }



    void test_local(init::QgodLocal::view::type QgodLocal, std::array<std::array<real, 9>, 9> solution) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          if (i == j) {
            test_relative(QgodLocal(i,j), 1-solution[j][i]);
          }
          else {
            test_relative(QgodLocal(i,j), -solution[j][i]);
          }
        }
      } 
    }

    void test_neighbor(init::QgodLocal::view::type QgodNeighbor, std::array<std::array<real, 9>, 9> solution) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          test_relative(QgodNeighbor(i,j), solution[j][i]);
        }
      }
    }

    void test_nan(init::QgodLocal::view::type QgodNeighbor) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          //CXXTEST TS_ASSERT_IS_NAN is broken
          TS_ASSERT(std::isnan(QgodNeighbor(i,j)));
        }
      }
    }

    void test_relative(real a, real b) {
      if (std::abs(a) < epsilon) {
        TS_ASSERT_DELTA(0, b, epsilon);
      } else if (std::abs(b) < epsilon) {
        TS_ASSERT_DELTA(0, a, epsilon);
      } else {
        TS_ASSERT_DELTA(a, b, std::abs(a)*1e2*epsilon);
        TS_ASSERT_DELTA(a, b, std::abs(b)*1e2*epsilon);
      }
    }

};
