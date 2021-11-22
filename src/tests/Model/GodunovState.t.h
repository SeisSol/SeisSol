#include <cmath>

#include <Equations/datastructures.hpp>
#include <Model/common.hpp>
#include <Equations/Setup.h>

#include "values.h"

template<typename T>
void test_matrix(seissol::init::QgodLocal::view::type QgodLocal, T solution, double epsilon) {
  //compute the Frobenius norms squared: ||QgodLocal - solution||^2 and ||solution||^2
  double frobDiffSquared= 0.0;
  double frobASquared = 0.0;
  for (unsigned int i = 0; i < solution[0].size(); i++) {
    for (unsigned int j = 0; j < solution.size(); j++) {
      const auto diff = (QgodLocal(i,j) - solution[j][i]);
      frobDiffSquared += diff*diff;
      frobASquared += solution[j][i]*solution[j][i];
    }
  }
  //Check whether ||QgodLocal - solution||^2 < epsilon^2 * ||solution||^2
  REQUIRE(std::abs(frobDiffSquared) < epsilon * frobASquared * epsilon * frobASquared);
}

void test_nan(seissol::init::QgodLocal::view::type QgodNeighbor) {
  for (unsigned int i = 0; i < 9; i++) {
    for (unsigned int j = 0; j < 9; j++) {
      REQUIRE(std::isnan(QgodNeighbor(i, j)));
    }
  }
}

TEST_CASE("Godunov state is correct") {
  constexpr real epsilon = 1e2 * std::numeric_limits<real>::epsilon();

  real localData[seissol::tensor::QgodLocal::size()];
  real neighborData[seissol::tensor::QgodLocal::size()];
  seissol::init::QgodLocal::view::type QgodLocal = seissol::init::QgodLocal::view::create(localData);
  seissol::init::QgodNeighbor::view::type QgodNeighbor = seissol::init::QgodNeighbor::view::create(neighborData);
  QgodLocal.setZero();
  QgodNeighbor.setZero();

  //test homogeneous material
#ifdef USE_ANISOTROPIC
  seissol::model::AnisotropicMaterial local(materialVal_1, 22);
      seissol::model::AnisotropicMaterial neighbor(materialVal_1, 22);
#elif defined USE_POROELASTIC
  seissol::model::PoroElasticMaterial local(materialVal_1, 10);
      seissol::model::PoroElasticMaterial neighbor(materialVal_1, 10);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  seissol::model::ViscoElasticMaterial local(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS*4);
      seissol::model::ViscoElasticMaterial neighbor(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS*4);
#else
  seissol::model::ElasticMaterial local(seissol::unit_test::materialVal_1, 3);
  seissol::model::ElasticMaterial neighbor(seissol::unit_test::materialVal_1, 3);
#endif

  seissol::model::getTransposedGodunovState( local,
                                             neighbor,
                                             FaceType::regular,
                                             QgodLocal,
                                             QgodNeighbor );
  test_matrix(QgodLocal, seissol::unit_test::solution_homogeneous_local, epsilon);
  test_matrix(QgodNeighbor, seissol::unit_test::solution_homogeneous_neighbor, epsilon);

  //test free surface
  seissol::model::getTransposedGodunovState( local,
                                             neighbor,
                                             FaceType::freeSurface,
                                             QgodLocal,
                                             QgodNeighbor );
  test_matrix(QgodLocal, seissol::unit_test::solution_boundary, epsilon);
  test_nan(QgodNeighbor);



  //test heterogeneous material
#ifdef USE_ANISOTROPIC
  neighbor = seissol::model::AnisotropicMaterial(materialVal_2, 22);
#elif defined USE_POROELASTIC
  neighbor = seissol::model::PoroElasticMaterial(materialVal_2, 10);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  neighbor = seissol::model::ViscoElasticMaterial(materialVal_2, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);
#else
  neighbor = seissol::model::ElasticMaterial(seissol::unit_test::materialVal_2, 3);
#endif

  seissol::model::getTransposedGodunovState( local,
                                             neighbor,
                                             FaceType::regular,
                                             QgodLocal,
                                             QgodNeighbor );
  test_matrix(QgodLocal, seissol::unit_test::solution_heterogeneous_local, epsilon);
  test_matrix(QgodNeighbor, seissol::unit_test::solution_heterogeneous_neighbor, epsilon);


}
