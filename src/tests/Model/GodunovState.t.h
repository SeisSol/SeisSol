// SPDX-FileCopyrightText: 2021 SeisSol Group
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
template <typename Cfg, typename T>
void test_matrix(typename init::QgodLocal<Cfg>::view::type& qgodLocal,
                 const T& solution,
                 double epsilon) {
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

template <typename Cfg>
inline void test_nan(typename init::QgodLocal<Cfg>::view::type& qgodNeighbor) {
  for (unsigned int i = 0; i < qgodNeighbor.shape(0); i++) {
    for (unsigned int j = 0; j < qgodNeighbor.shape(1); j++) {
      REQUIRE(std::isnan(qgodNeighbor(i, j)));
    }
  }
}

// TODO: add acoustic Godunov state
TEST_CASE_TEMPLATE_DEFINE("Godunov state is correct", Cfg, configId4) {
  using real = Real<Cfg>;
  using MaterialT = model::MaterialTT<Cfg>;

  if constexpr (std::is_same_v<MaterialT, model::AcousticMaterial>) {
    return;
  }

  constexpr real Epsilon = 1e2 * std::numeric_limits<real>::epsilon();

  real localData[tensor::QgodLocal<Cfg>::size()];
  real neighborData[tensor::QgodLocal<Cfg>::size()];
  typename init::QgodLocal<Cfg>::view::type qgodLocal =
      init::QgodLocal<Cfg>::view::create(localData);
  typename init::QgodNeighbor<Cfg>::view::type qgodNeighbor =
      init::QgodNeighbor<Cfg>::view::create(neighborData);
  qgodLocal.setZero();
  qgodNeighbor.setZero();

  // test homogeneous material
  const MaterialT local(SolutionData<MaterialT>::MaterialVal1);
  MaterialT neighbor(SolutionData<MaterialT>::MaterialVal1);

  seissol::model::getTransposedGodunovState<Cfg>(
      local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
  test_matrix<Cfg>(qgodLocal, SolutionData<MaterialT>::SolutionHomogeneousLocal, Epsilon);
  test_matrix<Cfg>(qgodNeighbor, SolutionData<MaterialT>::SolutionHomogeneousNeighbor, Epsilon);

  // test free surface
  seissol::model::getTransposedGodunovState<Cfg>(
      local, neighbor, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
  test_matrix<Cfg>(qgodLocal, SolutionData<MaterialT>::SolutionBoundary, Epsilon);
  test_nan<Cfg>(qgodNeighbor);

  // test heterogeneous material
  neighbor = MaterialT(SolutionData<MaterialT>::MaterialVal2);

  seissol::model::getTransposedGodunovState<Cfg>(
      local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
  test_matrix<Cfg>(qgodLocal, SolutionData<MaterialT>::SolutionHeterogeneousLocal, Epsilon);
  test_matrix<Cfg>(qgodNeighbor, SolutionData<MaterialT>::SolutionHeterogeneousNeighbor, Epsilon);
}
} // namespace seissol::unit_test
