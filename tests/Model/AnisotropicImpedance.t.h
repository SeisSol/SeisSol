// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_ANISOTROPICIMPEDANCE_T_H_
#define SEISSOL_TESTS_MODEL_ANISOTROPICIMPEDANCE_T_H_

// The whole file only makes sense for a build whose MaterialT is anisotropic.
#ifdef USE_ANISOTROPIC

#include <doctest.h>

#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Initializer/Model/DynamicRuptureImpedance.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <random>
#include <vector>

namespace seissol::unit_test {

namespace {

using seissol::initializer::model::checkFaultImpedance;
using seissol::initializer::model::computeAdmittance;
using seissol::initializer::model::computeFaultImpedance;
using seissol::initializer::model::DrMatrix;

struct FaultFrame {
  Eigen::Vector3d normal;
  Eigen::Vector3d tangent1;
  Eigen::Vector3d tangent2;
};

/// Right-handed orthonormal frame with the given (unnormalised) normal.
FaultFrame makeFaultFrame(const Eigen::Vector3d& rawNormal) {
  FaultFrame frame;
  frame.normal = rawNormal.normalized();
  const Eigen::Vector3d seed =
      (std::abs(frame.normal.x()) < 0.9) ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
  frame.tangent1 = frame.normal.cross(seed).normalized();
  frame.tangent2 = frame.normal.cross(frame.tangent1);
  return frame;
}

/// Rotates a material into the fault-local frame, exactly as
/// initializeDynamicRuptureMatrices does.
model::AnisotropicMaterial rotateToFault(const model::AnisotropicMaterial& material,
                                         const FaultFrame& frame) {
  const VrtxCoords normal{frame.normal.x(), frame.normal.y(), frame.normal.z()};
  const VrtxCoords tangent1{frame.tangent1.x(), frame.tangent1.y(), frame.tangent1.z()};
  const VrtxCoords tangent2{frame.tangent2.x(), frame.tangent2.y(), frame.tangent2.z()};
  std::array<double, 36> bond{};
  model::getBondMatrix(normal, tangent1, tangent2, bond);
  return model::getRotatedMaterialCoefficients(bond, material);
}

model::AnisotropicMaterial isotropicMaterial(double rho, double mu, double lambda) {
  return model::AnisotropicMaterial(model::ElasticMaterial(std::vector<double>{rho, mu, lambda}));
}

/// VTI with the symmetry axis along z, tilted by `tiltDegrees` towards x.
model::AnisotropicMaterial tiltedVti(double tiltDegrees) {
  model::AnisotropicMaterial vti;
  vti.rho = 2200.0;
  const double c11 = 3.10e10;
  const double c33 = 2.25e10;
  const double c13 = 1.06e10;
  const double c44 = 0.90e10;
  const double c66 = 1.15e10;
  vti.c11 = vti.c22 = c11;
  vti.c33 = c33;
  vti.c12 = c11 - 2 * c66;
  vti.c13 = vti.c23 = c13;
  vti.c44 = vti.c55 = c44;
  vti.c66 = c66;

  const double angle = tiltDegrees * M_PI / 180.0;
  // rotate the material about y, i.e. tilt the symmetry axis away from z
  const Eigen::Vector3d axis1(std::cos(angle), 0.0, -std::sin(angle));
  const Eigen::Vector3d axis2(0.0, 1.0, 0.0);
  const Eigen::Vector3d axis3(std::sin(angle), 0.0, std::cos(angle));
  const VrtxCoords a1{axis1.x(), axis1.y(), axis1.z()};
  const VrtxCoords a2{axis2.x(), axis2.y(), axis2.z()};
  const VrtxCoords a3{axis3.x(), axis3.y(), axis3.z()};
  std::array<double, 36> bond{};
  model::getBondMatrix(a1, a2, a3, bond);
  return model::getRotatedMaterialCoefficients(bond, vti);
}

double relError(const DrMatrix& value, const DrMatrix& reference) {
  return (value - reference).cwiseAbs().maxCoeff() / reference.cwiseAbs().maxCoeff();
}

} // namespace

// ---------------------------------------------------------------------------
// 1. Isotropic limit: an anisotropic material built from (lambda, mu) must
//    reproduce exactly the numbers the elastic "fast" branch of
//    initializeDynamicRuptureMatrices writes, for *any* fault orientation.
//    This is the regression test for the whole aniso DR impedance path.
// ---------------------------------------------------------------------------
TEST_CASE("Anisotropic DR impedance reduces to the elastic one" *
          doctest::test_suite("dynamicrupture")) {
  constexpr double Epsilon = 1e-12;

  const double rhoP = 2670.0;
  const double muP = 3.203e10;
  const double lambdaP = 3.204e10;
  const double rhoM = 2500.0;
  const double muM = 1.600e10;
  const double lambdaM = 2.000e10;

  const auto plus = isotropicMaterial(rhoP, muP, lambdaP);
  const auto minus = isotropicMaterial(rhoM, muM, lambdaM);

  const double zpP = rhoP * std::sqrt((lambdaP + 2 * muP) / rhoP);
  const double zsP = rhoP * std::sqrt(muP / rhoP);
  const double zpM = rhoM * std::sqrt((lambdaM + 2 * muM) / rhoM);
  const double zsM = rhoM * std::sqrt(muM / rhoM);
  const double etaP = zpP * zpM / (zpP + zpM);
  const double etaS = zsP * zsM / (zsP + zsM);

  const DrMatrix admittanceRefP = Eigen::Vector3d(1 / zpP, 1 / zsP, 1 / zsP).asDiagonal();
  const DrMatrix etaRef = Eigen::Vector3d(etaP, etaS, etaS).asDiagonal();
  // exactly the entries the elastic branch stores in tractionPlusMatrix
  const DrMatrix bRefP = Eigen::Vector3d(etaP / zpP, etaS / zsP, etaS / zsP).asDiagonal();
  const DrMatrix bRefM = Eigen::Vector3d(etaP / zpM, etaS / zsM, etaS / zsM).asDiagonal();

  // deterministic pseudo-random orientations plus the axis-aligned ones
  std::vector<Eigen::Vector3d> normals{Eigen::Vector3d::UnitX(),
                                       Eigen::Vector3d::UnitY(),
                                       Eigen::Vector3d::UnitZ(),
                                       Eigen::Vector3d(1, 1, 1),
                                       Eigen::Vector3d(0.3, -0.5, 0.81)};
  std::mt19937 rng(20260803);
  std::normal_distribution<double> normalDist;
  for (int i = 0; i < 64; ++i) {
    normals.emplace_back(normalDist(rng), normalDist(rng), normalDist(rng));
  }

  for (const auto& rawNormal : normals) {
    const auto frame = makeFaultFrame(rawNormal);
    const auto impedance =
        computeFaultImpedance(rotateToFault(plus, frame), rotateToFault(minus, frame));

    const auto violation = checkFaultImpedance(impedance);
    REQUIRE_MESSAGE(!violation.has_value(), violation.value_or(""));

    CHECK(relError(impedance.admittancePlus, admittanceRefP) < Epsilon);
    CHECK(relError(impedance.eta, etaRef) < Epsilon);
    CHECK(relError(impedance.bPlus, bRefP) < Epsilon);
    CHECK(relError(impedance.bMinus, bRefM) < Epsilon);
  }
}

// ---------------------------------------------------------------------------
// 2. Genuine anisotropy: invariants, and the normal/shear coupling that the
//    friction solvers have to pick up. The 0 deg / 90 deg cases are the guard
//    against a *missing* Bond rotation -- without it the coupling would not
//    depend on the fault orientation at all.
// ---------------------------------------------------------------------------
TEST_CASE("Anisotropic DR impedance has orientation dependent normal coupling" *
          doctest::test_suite("dynamicrupture")) {
  constexpr double Epsilon = 1e-12;

  const auto frame = makeFaultFrame(Eigen::Vector3d::UnitZ());

  const auto couplingRatio = [&](double tiltDegrees) {
    const auto material = rotateToFault(tiltedVti(tiltDegrees), frame);
    const auto impedance = computeFaultImpedance(material, material);
    const auto violation = checkFaultImpedance(impedance);
    REQUIRE_MESSAGE(!violation.has_value(), violation.value_or(""));
    // homogeneous fault: b+ must be exactly half the identity
    CHECK(relError(impedance.bPlus, DrMatrix(0.5 * DrMatrix::Identity())) < Epsilon);
    const double shear = std::max(std::abs(impedance.eta(1, 1)), std::abs(impedance.eta(2, 2)));
    return std::max(std::abs(impedance.eta(0, 1)), std::abs(impedance.eta(0, 2))) / shear;
  };

  // normal along / perpendicular to the symmetry axis: no normal-shear coupling
  CHECK(couplingRatio(0.0) < Epsilon);
  CHECK(couplingRatio(90.0) < Epsilon);

  // obliquely tilted symmetry axis: coupling of a few percent, peaking near 30 deg
  CHECK(couplingRatio(30.0) > 0.05);
  CHECK(couplingRatio(30.0) > couplingRatio(10.0));
  CHECK(couplingRatio(30.0) > couplingRatio(70.0));
}

// ---------------------------------------------------------------------------
// 3. Pins down the flat CSC layout of tractionPlusMatrix that
//    common::computeFrictionEnergy indexes directly (flat = 3 * col + row,
//    with row running over the *stored* rows {0, 3, 5}).
// ---------------------------------------------------------------------------
TEST_CASE("tractionPlusMatrix CSC layout matches the friction energy indexing" *
          doctest::test_suite("dynamicrupture")) {
  constexpr std::array<std::size_t, 3> StoredRows{0, 3, 5};
  constexpr std::size_t Rows = 3;

  REQUIRE(tensor::tractionPlusMatrix::size() == Rows * 3);

  alignas(Alignment) real data[tensor::tractionPlusMatrix::size()]{};
  auto view = init::tractionPlusMatrix::view::create(data);
  view.setZero();
  for (std::size_t col = 0; col < 3; ++col) {
    for (std::size_t row = 0; row < StoredRows.size(); ++row) {
      view(StoredRows[row], col) = static_cast<real>(10 * col + row);
    }
  }

  for (std::size_t col = 0; col < 3; ++col) {
    for (std::size_t row = 0; row < StoredRows.size(); ++row) {
      CHECK(data[Rows * col + row] == doctest::Approx(10.0 * col + row));
    }
  }
}

// ---------------------------------------------------------------------------
// 4. Reconstruction of the stress components outside the fault-normal Riemann problem, used by
//    the fault receiver output. Checked against the direct plane-wave stress, and against the
//    1 - 2 (cs/cp)^2 formula the output used before.
// ---------------------------------------------------------------------------
TEST_CASE("Anisotropic lateral stress reconstruction" * doctest::test_suite("dynamicrupture")) {
  using seissol::initializer::model::DrLateralMatrix;
  constexpr double Epsilon = 1e-12;

  // Voigt stiffness matrix of a material, for the direct check below
  const auto voigt = [](const model::AnisotropicMaterial& m) {
    Eigen::Matrix<double, 6, 6> c;
    c << m.c11, m.c12, m.c13, m.c14, m.c15, m.c16, m.c12, m.c22, m.c23, m.c24, m.c25, m.c26, m.c13,
        m.c23, m.c33, m.c34, m.c35, m.c36, m.c14, m.c24, m.c34, m.c44, m.c45, m.c46, m.c15, m.c25,
        m.c35, m.c45, m.c55, m.c56, m.c16, m.c26, m.c36, m.c46, m.c56, m.c66;
    return c;
  };

  SUBCASE("isotropic limit") {
    const double rho = 2670.0;
    const double mu = 3.203e10;
    const double lambda = 3.204e10;
    const auto frame = makeFaultFrame(Eigen::Vector3d(0.3, -0.5, 0.81));
    const auto lateral =
        seissol::initializer::model::impedance_detail::lateralStressFromChristoffel(
            rotateToFault(isotropicMaterial(rho, mu, lambda), frame));

    const double cs = std::sqrt(mu / rho);
    const double cp = std::sqrt((lambda + 2 * mu) / rho);
    Eigen::Matrix3d reference = Eigen::Matrix3d::Zero();
    reference(0, 0) = 1.0 - 2.0 * (cs * cs) / (cp * cp);
    reference(1, 0) = reference(0, 0);

    CHECK((lateral - reference).cwiseAbs().maxCoeff() < Epsilon);
  }

  SUBCASE("reproduces the plane wave stress of a tilted VTI") {
    const auto frame = makeFaultFrame(Eigen::Vector3d::UnitZ());
    const auto material = rotateToFault(tiltedVti(35.0), frame);
    const auto lateral =
        seissol::initializer::model::impedance_detail::lateralStressFromChristoffel(material);
    const auto stiffness = voigt(material);

    // every polarisation of a wave travelling along the fault normal has only the strain rates
    // (eps_1, eps_6, eps_5); its full stress must be recovered from its traction alone
    Eigen::Matrix3d christoffel;
    christoffel << material.c11, material.c16, material.c15, material.c16, material.c66,
        material.c56, material.c15, material.c56, material.c55;
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(christoffel);

    for (int mode = 0; mode < 3; ++mode) {
      const Eigen::Vector3d polarisation = solver.eigenvectors().col(mode);
      Eigen::Matrix<double, 6, 1> strain = Eigen::Matrix<double, 6, 1>::Zero();
      strain(0) = polarisation(0);
      strain(5) = polarisation(1);
      strain(4) = polarisation(2);
      const Eigen::Matrix<double, 6, 1> stress = stiffness * strain;

      const Eigen::Vector3d traction(stress(0), stress(5), stress(4));
      const Eigen::Vector3d expected(stress(1), stress(2), stress(3));

      CHECK((lateral * traction - expected).cwiseAbs().maxCoeff() <
            Epsilon * expected.cwiseAbs().maxCoeff());
    }

    // the shear-to-lateral coupling is absent from the isotropic formula, and it is not small
    CHECK(lateral.rightCols<2>().cwiseAbs().maxCoeff() > 0.05 * lateral.cwiseAbs().maxCoeff());
  }

  SUBCASE("closed form agrees with the eigendecomposition route") {
    // The poroelastic case has no closed form and goes through
    // admittanceFromEigendecomposition instead. Checking the two against each other here -- on a
    // material where the eigensolver is well conditioned -- validates the row/column bookkeeping
    // of that route, which is all poroelasticity relies on.
    const auto frame = makeFaultFrame(Eigen::Vector3d(0.3, -0.5, 0.81));
    const auto material = rotateToFault(tiltedVti(35.0), frame);

    DrLateralMatrix fromEigen = DrLateralMatrix::Zero();
    seissol::initializer::model::impedance_detail::admittanceFromEigendecomposition(material,
                                                                                    &fromEigen);
    const auto fromClosedForm =
        seissol::initializer::model::impedance_detail::lateralStressFromChristoffel(material);

    CHECK((fromEigen - fromClosedForm).cwiseAbs().maxCoeff() <
          1e-8 * fromClosedForm.cwiseAbs().maxCoeff());
  }
}

} // namespace seissol::unit_test

#endif // USE_ANISOTROPIC

#endif // SEISSOL_TESTS_MODEL_ANISOTROPICIMPEDANCE_T_H_
