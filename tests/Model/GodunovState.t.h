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
#include <vector>

namespace seissol::unit_test {

// silence double == real cast
inline double castReal(real value) {
  if constexpr (std::is_same_v<double, real>) {
    return value;
  } else {
    return static_cast<double>(value);
  }
}

/**
 * Checks ||A - B||_F <= epsilon * ||B||_F, i.e. a relative comparison in the Frobenius norm.
 *
 * Both arguments are the squared norms, so the check is ||A - B||_F^2 <= epsilon^2 * ||B||_F^2.
 * Getting this wrong is easy and expensive: comparing against (epsilon * ||B||_F^2)^2 -- as this
 * file used to -- turns the bound into ||A - B||_F <= epsilon * ||B||_F^2, which for the Godunov
 * matrices (entries of order 1e6, ||B||_F^2 ~ 1e14) permits an absolute deviation of ~2.8 and makes
 * the assertion vacuous.
 */
inline void checkRelative(double frobDiffSquared, double frobASquared, double epsilon) {
  CHECK(frobDiffSquared <= epsilon * epsilon * frobASquared);
}

template <typename T>
void testMatrix(init::QgodLocal::view::type& qgod, const T& solution, double epsilon) {
  double frobDiffSquared = 0.0;
  double frobASquared = 0.0;
  for (std::size_t i = 0; i < solution[0].size(); i++) {
    for (std::size_t j = 0; j < solution.size(); j++) {
      const auto diff = castReal(qgod(i, j)) - solution[j][i];
      frobDiffSquared += diff * diff;
      frobASquared += solution[j][i] * solution[j][i];
    }
  }
  checkRelative(frobDiffSquared, frobASquared, epsilon);
}

/**
 * QgodLocal + QgodNeighbor == I.
 *
 * This is an invariant of the Godunov decomposition, not of any particular implementation: the two
 * selectors partition the full set of characteristic modes, so the two projectors are complementary
 * by construction. The current setups derive QgodLocal as I - QgodNeighbor and therefore satisfy it
 * to the last bit; the check guards against a regression to two independent solves, which would
 * leave a residual of order cond(matR) * eps (~1e-9 for realistic material contrasts), and against
 * a mode partition that leaves the null space out of both matrices (which is what the poroelastic
 * setup used to do).
 */
inline void testConsistency(init::QgodLocal::view::type& qgodLocal,
                            init::QgodNeighbor::view::type& qgodNeighbor,
                            double epsilon) {
  double diffSquared = 0.0;
  double frobASquared = 0.0;
  for (std::size_t i = 0; i < qgodNeighbor.shape(0); i++) {
    for (std::size_t j = 0; j < qgodNeighbor.shape(1); j++) {
      const auto sol = (i == j) ? 1.0 : 0.0;
      const auto diff = (castReal(qgodLocal(i, j)) + castReal(qgodNeighbor(i, j))) - sol;
      diffSquared += diff * diff;
      frobASquared += sol * sol;
    }
  }
  checkRelative(diffSquared, frobASquared, epsilon);
}

/**
 * QgodNeighbor is a projector, i.e. P^2 == P.
 *
 * Unlike the consistency check this cannot be satisfied by construction -- it is a genuine
 * numerical statement about the quality of the eigenbasis that matR was built from. It is the
 * diagnostic that detects a (near-)degenerate eigenvector basis: with the eigenvectors that a
 * general eigensolver returns for the doubly degenerate poroelastic shear eigenvalue this
 * residual sits at ~7e-12 in double, far above the ~1e-16 that a well-conditioned basis achieves.
 *
 * ||P^2 - P||_F is compared against ||P||_F (not ||P||_F^2): squaring a matrix with entries of
 * order 1e6 that cancel back down to order 1e6 amplifies the input rounding by ||P||_F, so this is
 * the scaling at which the residual is precision-independent (~1e-16 in double, ~1e-8 in single).
 */
inline void testProjector(init::QgodNeighbor::view::type& qgodNeighbor, double epsilon) {
  const std::size_t rows = qgodNeighbor.shape(0);
  const std::size_t cols = qgodNeighbor.shape(1);

  std::vector<double> matP(rows * cols);
  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = 0; j < cols; j++) {
      matP[i * cols + j] = castReal(qgodNeighbor(i, j));
    }
  }

  double frobDiffSquared = 0.0;
  double frobASquared = 0.0;
  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = 0; j < cols; j++) {
      double squared = 0.0;
      for (std::size_t k = 0; k < cols; k++) {
        squared += matP[i * cols + k] * matP[k * cols + j];
      }
      const auto diff = squared - matP[i * cols + j];
      frobDiffSquared += diff * diff;
      frobASquared += matP[i * cols + j] * matP[i * cols + j];
    }
  }
  checkRelative(frobDiffSquared, frobASquared, epsilon);
}

inline void testNAN(init::QgodNeighbor::view::type& qgodNeighbor) {
  for (std::size_t i = 0; i < qgodNeighbor.shape(0); i++) {
    for (std::size_t j = 0; j < qgodNeighbor.shape(1); j++) {
      CHECK(std::isnan(qgodNeighbor(i, j)));
    }
  }
}

TEST_CASE("Godunov state is correct" * doctest::test_suite("model")) {
  // Tolerance for comparing against the stored reference matrices. This is a property of the
  // equation set (see SolutionData::MatrixEpsilon): the elastic family assembles its eigenbasis in
  // closed form and hits machine precision, whereas the poroelastic one goes through a numerical
  // eigendecomposition and cannot.
  constexpr double MatrixEpsilon = SolutionData<model::MaterialT>::MatrixEpsilon;

  // The structural invariants below hold independently of how well conditioned matR is, so they are
  // checked at machine precision for every equation set.
  constexpr double StructuralEpsilon = 1e2 * std::numeric_limits<real>::epsilon();

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
    testMatrix(qgodLocal, SolutionData<model::MaterialT>::SolutionHomogeneousLocal, MatrixEpsilon);
    testMatrix(
        qgodNeighbor, SolutionData<model::MaterialT>::SolutionHomogeneousNeighbor, MatrixEpsilon);
    testConsistency(qgodLocal, qgodNeighbor, StructuralEpsilon);
    testProjector(qgodNeighbor, StructuralEpsilon);
  }

  SUBCASE("Free surface material") {
    // material 1 vs 1
    const model::MaterialT local(SolutionData<model::MaterialT>::MaterialVal1);
    const model::MaterialT neighbor(SolutionData<model::MaterialT>::MaterialVal1);

    model::getTransposedGodunovState(
        local, neighbor, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
    testMatrix(qgodLocal, SolutionData<model::MaterialT>::SolutionBoundary, MatrixEpsilon);
    testNAN(qgodNeighbor);
  }

  SUBCASE("Heterogenous material") {
    // material 1 vs 2
    const model::MaterialT local(SolutionData<model::MaterialT>::MaterialVal1);
    const model::MaterialT neighbor(SolutionData<model::MaterialT>::MaterialVal2);

    model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
    testMatrix(
        qgodLocal, SolutionData<model::MaterialT>::SolutionHeterogeneousLocal, MatrixEpsilon);
    testMatrix(
        qgodNeighbor, SolutionData<model::MaterialT>::SolutionHeterogeneousNeighbor, MatrixEpsilon);
    testConsistency(qgodLocal, qgodNeighbor, StructuralEpsilon);
    testProjector(qgodNeighbor, StructuralEpsilon);
  }
}
} // namespace seissol::unit_test
