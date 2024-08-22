#include <cmath>

#include "Equations/Setup.h"
#include "Equations/datastructures.hpp"
#include "Model/common.hpp"

#include "values.h"

namespace seissol::unit_test {
template <typename T>
void test_matrix(init::QgodLocal::view::type qgodLocal, T solution, double epsilon) {
  // compute the Frobenius norms squared: ||QgodLocal - solution||^2 and ||solution||^2
  double frobDiffSquared = 0.0;
  double frobASquared = 0.0;
  for (unsigned int i = 0; i < solution[0].size(); i++) {
    for (unsigned int j = 0; j < solution.size(); j++) {
      const auto diff = (qgodLocal(i, j) - solution[j][i]);
      frobDiffSquared += diff * diff;
      frobASquared += solution[j][i] * solution[j][i];
    }
  }
  // Check whether ||QgodLocal - solution||^2 < epsilon^2 * ||solution||^2
  REQUIRE(std::abs(frobDiffSquared) < epsilon * frobASquared * epsilon * frobASquared);
}

inline void test_nan(init::QgodLocal::view::type qgodNeighbor) {
  for (unsigned int i = 0; i < 9; i++) {
    for (unsigned int j = 0; j < 9; j++) {
      REQUIRE(std::isnan(qgodNeighbor(i, j)));
    }
  }
}

TEST_CASE("Godunov state is correct") {
  constexpr real Epsilon = 1e2 * std::numeric_limits<real>::epsilon();

  real localData[tensor::QgodLocal::size()];
  real neighborData[tensor::QgodLocal::size()];
  init::QgodLocal::view::type qgodLocal = init::QgodLocal::view::create(localData);
  init::QgodNeighbor::view::type qgodNeighbor = init::QgodNeighbor::view::create(neighborData);
  qgodLocal.setZero();
  qgodNeighbor.setZero();

  // test homogeneous material
#ifdef USE_ANISOTROPIC
  model::AnisotropicMaterial local(materialVal_1, 22);
  model::AnisotropicMaterial neighbor(materialVal_1, 22);
#elif defined USE_POROELASTIC
  model::PoroElasticMaterial local(materialVal_1, 10);
  model::PoroElasticMaterial neighbor(materialVal_1, 10);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  model::ViscoElasticMaterial local(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);
  model::ViscoElasticMaterial neighbor(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);
#else
  const model::ElasticMaterial local(materialVal1, 3);
  model::ElasticMaterial neighbor(materialVal1, 3);
#endif

  model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
  test_matrix(qgodLocal, solutionHomogeneousLocal, Epsilon);
  test_matrix(qgodNeighbor, solutionHomogeneousNeighbor, Epsilon);

  // test free surface
  model::getTransposedGodunovState(local, neighbor, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
  test_matrix(qgodLocal, unit_test::solutionBoundary, Epsilon);
  test_nan(qgodNeighbor);

  // test heterogeneous material
#ifdef USE_ANISOTROPIC
  neighbor = model::AnisotropicMaterial(materialVal_2, 22);
#elif defined USE_POROELASTIC
  neighbor = model::PoroElasticMaterial(materialVal_2, 10);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  neighbor = model::ViscoElasticMaterial(materialVal_2, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);
#else
  neighbor = model::ElasticMaterial(materialVal2, 3);
#endif

  model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
  test_matrix(qgodLocal, solutionHeterogeneousLocal, Epsilon);
  test_matrix(qgodNeighbor, solutionHeterogeneousNeighbor, Epsilon);
}
} // namespace seissol::unit_test
