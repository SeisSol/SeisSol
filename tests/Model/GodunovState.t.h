// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "Equations/acoustic/Model/Datastructures.h"
#include "Model/Common.h"
#include "TestHelper.h"
#include "Values.h"

#include <cmath>

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
  CHECK(std::abs(frobDiffSquared) < epsilon * frobASquared * epsilon * frobASquared);
}

inline void test_consistency(init::QgodLocal::view::type& qgodLocal,
                             init::QgodLocal::view::type& qgodNeighbor,
                             double epsilon) {
  // compute the Frobenius norms squared: ||QgodLocal - solution||^2 and ||solution||^2
  double diffSquared = 0.0;
  double frobASquared = 0.0;
  for (unsigned int i = 0; i < qgodNeighbor.shape(0); i++) {
    for (unsigned int j = 0; j < qgodNeighbor.shape(1); j++) {
      const auto sol = i == j ? 1.0 : 0.0;
      const auto diff = (qgodLocal(i, j) + qgodNeighbor(i, j)) - sol;
      diffSquared += diff * diff;
      frobASquared += sol;
    }
  }
  // Check whether ||QgodLocal - solution||^2 < epsilon^2 * ||solution||^2
  CHECK(std::abs(diffSquared) < epsilon * frobASquared * epsilon * frobASquared);
}

inline void test_nan(init::QgodLocal::view::type& qgodNeighbor) {
  for (unsigned int i = 0; i < qgodNeighbor.shape(0); i++) {
    for (unsigned int j = 0; j < qgodNeighbor.shape(1); j++) {
      CHECK(std::isnan(qgodNeighbor(i, j)));
    }
  }
}

TEST_CASE("Godunov state is correct" * doctest::test_suite("model")) {
  constexpr real Epsilon = 1e2 * std::numeric_limits<real>::epsilon();

  real localData[tensor::QgodLocal::size()]{};
  real neighborData[tensor::QgodNeighbor::size()]{};
  init::QgodLocal::view::type qgodLocal = init::QgodLocal::view::create(localData);
  init::QgodNeighbor::view::type qgodNeighbor = init::QgodNeighbor::view::create(neighborData);
  qgodLocal.setZero();
  qgodNeighbor.setZero();

  SUBCASE("Homogenous material") {
    // material 1 vs 1
    const model::MaterialT local(SolutionData<model::MaterialT>::MaterialVal1);
    const model::MaterialT neighbor(SolutionData<model::MaterialT>::MaterialVal1);

    model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
    test_matrix(qgodLocal, SolutionData<model::MaterialT>::SolutionHomogeneousLocal, Epsilon);
    test_matrix(qgodNeighbor, SolutionData<model::MaterialT>::SolutionHomogeneousNeighbor, Epsilon);
    test_consistency(qgodLocal, qgodNeighbor, Epsilon);
  }

  SUBCASE("Free surface material") {
    // material 1 vs 1
    const model::MaterialT local(SolutionData<model::MaterialT>::MaterialVal1);
    const model::MaterialT neighbor(SolutionData<model::MaterialT>::MaterialVal1);

    model::getTransposedGodunovState(
        local, neighbor, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
    test_matrix(qgodLocal, SolutionData<model::MaterialT>::SolutionBoundary, Epsilon);
    test_nan(qgodNeighbor);
  }

  SUBCASE("Heterogenous material") {
    // material 1 vs 2
    const model::MaterialT local(SolutionData<model::MaterialT>::MaterialVal1);
    const model::MaterialT neighbor(SolutionData<model::MaterialT>::MaterialVal2);

    model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
    test_matrix(qgodLocal, SolutionData<model::MaterialT>::SolutionHeterogeneousLocal, Epsilon);
    test_matrix(
        qgodNeighbor, SolutionData<model::MaterialT>::SolutionHeterogeneousNeighbor, Epsilon);
    test_consistency(qgodLocal, qgodNeighbor, Epsilon);
  }
}
} // namespace seissol::unit_test
