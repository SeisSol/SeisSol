// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <Equations/acoustic/Model/Datastructures.h>
#include <cmath>
#include <tests/TestHelper.h>

#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "Model/Common.h"

#include "Values.h"

namespace seissol::unit_test {
template <typename T>
void test_matrix(init::QgodLocal::view::type& qgodLocal, const T& solution, double epsilon) {
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

inline void test_nan(init::QgodLocal::view::type& qgodNeighbor) {
  for (unsigned int i = 0; i < qgodNeighbor.shape(0); i++) {
    for (unsigned int j = 0; j < qgodNeighbor.shape(1); j++) {
      REQUIRE(std::isnan(qgodNeighbor(i, j)));
    }
  }
}

// TODO: add acoustic Godunov state
TEST_CASE("Godunov state is correct" *
          doctest::skip(std::is_same_v<model::MaterialT, model::AcousticMaterial>)) {
  constexpr real Epsilon = 1e2 * std::numeric_limits<real>::epsilon();

  real localData[tensor::QgodLocal::size()];
  real neighborData[tensor::QgodLocal::size()];
  init::QgodLocal::view::type qgodLocal = init::QgodLocal::view::create(localData);
  init::QgodNeighbor::view::type qgodNeighbor = init::QgodNeighbor::view::create(neighborData);
  qgodLocal.setZero();
  qgodNeighbor.setZero();

  // test homogeneous material
  const model::MaterialT local(SolutionData<model::MaterialT>::MaterialVal1);
  model::MaterialT neighbor(SolutionData<model::MaterialT>::MaterialVal1);

  model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
  test_matrix(qgodLocal, SolutionData<model::MaterialT>::SolutionHomogeneousLocal, Epsilon);
  test_matrix(qgodNeighbor, SolutionData<model::MaterialT>::SolutionHomogeneousNeighbor, Epsilon);

  // test free surface
  model::getTransposedGodunovState(local, neighbor, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
  test_matrix(qgodLocal, SolutionData<model::MaterialT>::SolutionBoundary, Epsilon);
  test_nan(qgodNeighbor);

  // test heterogeneous material
  neighbor = model::MaterialT(SolutionData<model::MaterialT>::MaterialVal2);

  model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
  test_matrix(qgodLocal, SolutionData<model::MaterialT>::SolutionHeterogeneousLocal, Epsilon);
  test_matrix(qgodNeighbor, SolutionData<model::MaterialT>::SolutionHeterogeneousNeighbor, Epsilon);
}
} // namespace seissol::unit_test
