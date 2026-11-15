// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "Equations/acoustic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Model/Common.h"
#include "Numerical/Eigenvalues.h"
#include "TestConfigs.h"
#include "TestHelper.h"
#include "Values.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>
#include <yateto.h>

namespace seissol::unit_test {

// silence double == real cast
template <typename RealT>
double castReal(RealT value) {
  if constexpr (std::is_same_v<double, RealT>) {
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

template <typename ViewT, typename T>
void testMatrix(ViewT& qgod, const T& solution, double epsilon) {
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
 * QgodLocal + QgodNeighbor == I, on the `quantities` the Godunov state acts on.
 *
 * This is an invariant of the Godunov decomposition, not of any particular implementation: the two
 * selectors partition the full set of characteristic modes, so the two projectors are complementary
 * by construction. The current setups derive QgodLocal as I - QgodNeighbor and therefore satisfy it
 * to the last bit; the check guards against a regression to two independent solves, which would
 * leave a residual of order cond(matR) * eps (~1e-9 for realistic material contrasts), and against
 * a mode partition that leaves the null space out of both matrices (which is what the poroelastic
 * setup used to do).
 */
template <typename LocalViewT, typename NeighborViewT>
void testConsistency(LocalViewT& qgodLocal,
                     NeighborViewT& qgodNeighbor,
                     std::size_t quantities,
                     double epsilon) {
  double diffSquared = 0.0;
  double frobASquared = 0.0;
  for (std::size_t i = 0; i < quantities; i++) {
    for (std::size_t j = 0; j < quantities; j++) {
      const auto sol = (i == j) ? 1.0 : 0.0;
      const auto diff = (castReal(qgodLocal(i, j)) + castReal(qgodNeighbor(i, j))) - sol;
      diffSquared += diff * diff;
      frobASquared += sol * sol;
    }
  }
  checkRelative(diffSquared, frobASquared, epsilon);
}

/**
 * QgodNeighbor is a projector, i.e. P^2 == P, on the `quantities` the Godunov state acts on.
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
template <typename ViewT>
void testProjector(ViewT& qgodNeighbor, std::size_t quantities, double epsilon) {
  const std::size_t rows = quantities;
  const std::size_t cols = quantities;

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

template <typename ViewT>
void testNAN(ViewT& qgodNeighbor, std::size_t quantities) {
  for (std::size_t i = 0; i < quantities; i++) {
    for (std::size_t j = 0; j < quantities; j++) {
      CHECK(std::isnan(qgodNeighbor(i, j)));
    }
  }
}

TEST_CASE_TEMPLATE("Godunov state is correct" * doctest::test_suite("model"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using MaterialT = model::MaterialOf<Cfg>;
  using real = Real<Cfg>;

  // Tolerance for comparing against the stored reference matrices. This is a property of the
  // equation set (see SolutionData::MatrixEpsilon): the elastic family assembles its eigenbasis in
  // closed form and hits machine precision, whereas the poroelastic one goes through a numerical
  // eigendecomposition and cannot.
  constexpr double MatrixEpsilon = SolutionData<MaterialT>::template MatrixEpsilon<real>;

  // The structural invariants below hold independently of how well conditioned matR is, so they are
  // checked at machine precision for every equation set.
  constexpr double StructuralEpsilon = 1e2 * std::numeric_limits<real>::epsilon();

  // The Godunov state acts on the quantities of the Riemann problem, i.e. the elastic ones. A
  // solver that keeps the anelastic quantities in the same tensor (fused mechanisms) only stores
  // the rows of these, so the checks stay within them.
  constexpr std::size_t Quantities = MaterialT::NumElasticQuantities;

  real localData[tensor::QgodLocal<Cfg>::size()]{};
  real neighborData[tensor::QgodNeighbor<Cfg>::size()]{};
  auto qgodLocal = init::QgodLocal<Cfg>::view::create(localData);
  auto qgodNeighbor = init::QgodNeighbor<Cfg>::view::create(neighborData);
  qgodLocal.setZero();
  qgodNeighbor.setZero();

  SUBCASE("Homogenous material") {
    // material 1 vs 1
    const MaterialT local(SolutionData<MaterialT>::MaterialVal1);
    const MaterialT neighbor(SolutionData<MaterialT>::MaterialVal1);

    model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
    testMatrix(qgodLocal, SolutionData<MaterialT>::SolutionHomogeneousLocal, MatrixEpsilon);
    testMatrix(qgodNeighbor, SolutionData<MaterialT>::SolutionHomogeneousNeighbor, MatrixEpsilon);
    testConsistency(qgodLocal, qgodNeighbor, Quantities, StructuralEpsilon);
    testProjector(qgodNeighbor, Quantities, StructuralEpsilon);
  }

  SUBCASE("Free surface material") {
    // material 1 vs 1
    const MaterialT local(SolutionData<MaterialT>::MaterialVal1);
    const MaterialT neighbor(SolutionData<MaterialT>::MaterialVal1);

    model::getTransposedGodunovState(
        local, neighbor, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
    testMatrix(qgodLocal, SolutionData<MaterialT>::SolutionBoundary, MatrixEpsilon);
    testNAN(qgodNeighbor, Quantities);
  }

  SUBCASE("Heterogenous material") {
    // material 1 vs 2
    const MaterialT local(SolutionData<MaterialT>::MaterialVal1);
    const MaterialT neighbor(SolutionData<MaterialT>::MaterialVal2);

    model::getTransposedGodunovState(local, neighbor, FaceType::Regular, qgodLocal, qgodNeighbor);
    testMatrix(qgodLocal, SolutionData<MaterialT>::SolutionHeterogeneousLocal, MatrixEpsilon);
    testMatrix(qgodNeighbor, SolutionData<MaterialT>::SolutionHeterogeneousNeighbor, MatrixEpsilon);
    testConsistency(qgodLocal, qgodNeighbor, Quantities, StructuralEpsilon);
    testProjector(qgodNeighbor, Quantities, StructuralEpsilon);
  }
}
TEST_CASE_TEMPLATE("Between an anisotropic and an isotropic solid, the Godunov state is the "
                   "isotropic one where the anisotropic solid is isotropic" *
                       doctest::test_suite("model"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using MaterialT = model::MaterialOf<Cfg>;
  if constexpr (model::SameFamily<MaterialT, model::ElasticMaterial>) {
    using real = Real<Cfg>;
    constexpr double Epsilon = SolutionData<MaterialT>::template MatrixEpsilon<real>;
    constexpr std::size_t Quantities = MaterialT::NumElasticQuantities;

    // the face with the normal along x, whose coordinates are the global ones
    std::array<double, 36> rotation{};
    model::getBondMatrix({1, 0, 0}, {0, 1, 0}, {0, 0, 1}, rotation);

    const model::ElasticMaterial elastic(SolutionData<model::ElasticMaterial>::MaterialVal1);
    const model::ElasticMaterial neighbor(SolutionData<model::ElasticMaterial>::MaterialVal2);

    real referenceLocalData[tensor::QgodLocal<Cfg>::size()]{};
    real referenceNeighborData[tensor::QgodNeighbor<Cfg>::size()]{};
    real localData[tensor::QgodLocal<Cfg>::size()]{};
    real neighborData[tensor::QgodNeighbor<Cfg>::size()]{};
    auto referenceLocal = init::QgodLocal<Cfg>::view::create(referenceLocalData);
    auto referenceNeighbor = init::QgodNeighbor<Cfg>::view::create(referenceNeighborData);
    auto qgodLocal = init::QgodLocal<Cfg>::view::create(localData);
    auto qgodNeighbor = init::QgodNeighbor<Cfg>::view::create(neighborData);

    const auto check = [&](auto& computed, auto& reference) {
      double diffSquared = 0.0;
      double referenceSquared = 0.0;
      for (std::size_t i = 0; i < Quantities; ++i) {
        for (std::size_t j = 0; j < Quantities; ++j) {
          const auto diff = castReal(computed(i, j)) - castReal(reference(i, j));
          diffSquared += diff * diff;
          referenceSquared += castReal(reference(i, j)) * castReal(reference(i, j));
        }
      }
      checkRelative(diffSquared, referenceSquared, Epsilon);
    };

    if constexpr (std::is_same_v<MaterialT, model::AnisotropicMaterial>) {
      // an anisotropic cell with an isotropic neighbor
      model::getTransposedGodunovState(
          elastic, neighbor, FaceType::Regular, referenceLocal, referenceNeighbor);
      model::getFaceGodunovState(rotation,
                                 model::AnisotropicMaterial(elastic),
                                 neighbor,
                                 FaceType::Regular,
                                 qgodLocal,
                                 qgodNeighbor);
    } else {
      // a cell of an isotropic solid with an anisotropic neighbor; the Riemann problem of a
      // viscoelastic solid is posed with its elastic part
      const MaterialT local(SolutionData<MaterialT>::MaterialVal1);
      model::getTransposedGodunovState(static_cast<const model::ElasticMaterial&>(local),
                                       neighbor,
                                       FaceType::Regular,
                                       referenceLocal,
                                       referenceNeighbor);
      model::getFaceGodunovState(rotation,
                                 local,
                                 model::AnisotropicMaterial(neighbor),
                                 FaceType::Regular,
                                 qgodLocal,
                                 qgodNeighbor);
    }
    check(qgodLocal, referenceLocal);
    check(qgodNeighbor, referenceNeighbor);
  }
}

/**
 * The Godunov states of a poroelastic solid hold its waves, whatever basis the eigensolver returns
 * for them. On a face within the solid, the projector P onto the outgoing waves (the transpose of
 * QgodNeighbor) commutes with the Jacobian A along the normal; at a free surface, the state maps
 * the tractions of every outgoing wave, a column of P, to minus its velocities. Both are measured
 * in the balanced basis of A (see eigenvalues::balance), in which the stresses and the velocities
 * of a wave compare. Besides the materials of the other tests, the materials are four, drawn at
 * random, for which Eigen 3.4.0 returns the eigenvectors of the shear waves parallel, or with
 * complex factors that leave their real parts close to zero.
 */
TEST_CASE_TEMPLATE("The Godunov states of a poroelastic solid hold its waves" *
                       doctest::test_suite("model"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using MaterialT = model::MaterialOf<Cfg>;
  if constexpr (std::is_same_v<MaterialT, model::PoroElasticMaterial>) {
    using real = Real<Cfg>;
    constexpr std::size_t Quantities = MaterialT::NumQuantities;
    constexpr double Epsilon = std::max(1e-13, 1e2 * std::numeric_limits<real>::epsilon());
    using Matrix = Eigen::Matrix<double, Quantities, Quantities>;

    const std::array<int, 4> tractions = {0, 3, 5, 9};
    const std::array<int, 6> velocities = {6, 7, 8, 10, 11, 12};
    const std::vector<std::vector<double>> materials = {SolutionData<MaterialT>::MaterialVal1,
                                                        SolutionData<MaterialT>::MaterialVal2,
                                                        {44564636678.115028,
                                                         2890.0262629452109,
                                                         21791374128.141182,
                                                         13692636429.049688,
                                                         0.073859362721722624,
                                                         1.8347793308397006e-15,
                                                         2.8557363252934556,
                                                         2477306928.3512659,
                                                         1005.7302325633887,
                                                         0.001},
                                                        {23779612930.962402,
                                                         2391.5765046602664,
                                                         8116829825.9581985,
                                                         17511912363.283329,
                                                         0.063505100561300445,
                                                         1.5665912677342567e-13,
                                                         2.3917663307696793,
                                                         2219620186.7326226,
                                                         1040.5566790854302,
                                                         0.001},
                                                        {35309325206.656868,
                                                         2554.3280264548694,
                                                         20172166757.33548,
                                                         2443508260.8886094,
                                                         0.34109863490187897,
                                                         4.2115433437873606e-15,
                                                         1.6533820054665773,
                                                         2242569169.0628324,
                                                         1048.9297326676815,
                                                         0.001},
                                                        {50339484502.353134,
                                                         2078.7721615636492,
                                                         20256335086.400734,
                                                         14144301454.696543,
                                                         0.33633906593111418,
                                                         2.5463020317398182e-16,
                                                         2.4666003072817011,
                                                         2317561307.7610407,
                                                         1028.5550706784359,
                                                         0.001}};

    for (const auto& values : materials) {
      const MaterialT material(values);

      std::array<double, Quantities * Quantities> transposedData{};
      auto transposed =
          yateto::DenseTensorView<2, double>(transposedData.data(), {Quantities, Quantities});
      model::MaterialSetup<MaterialT>::getTransposedCoefficientMatrix(material, 0, transposed);
      Matrix balanced;
      for (std::size_t i = 0; i < Quantities; ++i) {
        for (std::size_t j = 0; j < Quantities; ++j) {
          balanced(i, j) = transposed(j, i);
        }
      }
      const auto scale = eigenvalues::balance(balanced);

      real localData[tensor::QgodLocal<Cfg>::size()]{};
      real neighborData[tensor::QgodNeighbor<Cfg>::size()]{};
      auto qgodLocal = init::QgodLocal<Cfg>::view::create(localData);
      auto qgodNeighbor = init::QgodNeighbor<Cfg>::view::create(neighborData);
      model::getTransposedGodunovState(
          material, material, FaceType::Regular, qgodLocal, qgodNeighbor);
      Matrix projector;
      for (std::size_t i = 0; i < Quantities; ++i) {
        for (std::size_t j = 0; j < Quantities; ++j) {
          projector(i, j) = castReal(qgodNeighbor(j, i)) * scale(j) / scale(i);
        }
      }
      const double commutator = (balanced * projector - projector * balanced).cwiseAbs().maxCoeff();
      CHECK(commutator <=
            Epsilon * balanced.cwiseAbs().maxCoeff() * projector.cwiseAbs().maxCoeff());

      model::getTransposedGodunovState(
          material, material, FaceType::FreeSurface, qgodLocal, qgodNeighbor);
      for (std::size_t column = 0; column < Quantities; ++column) {
        double residual = 0.0;
        for (const auto velocity : velocities) {
          double sum = projector(velocity, column);
          for (const auto traction : tractions) {
            sum += castReal(qgodLocal(traction, velocity)) * scale(traction) / scale(velocity) *
                   projector(traction, column);
          }
          residual = std::max(residual, std::abs(sum));
        }
        CHECK(residual <= Epsilon * std::max(projector.col(column).cwiseAbs().maxCoeff(), 1.0));
      }
    }
  }
}

} // namespace seissol::unit_test
