// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/Eigenvalues.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <utility>

namespace seissol::unit_test {

template <std::size_t Dim>
void testResidual(std::array<std::complex<double>, Dim * Dim>& m,
                  seissol::eigenvalues::Eigenpair<std::complex<double>, Dim>& eigenpair) {
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 1e3;

  // compute M*R
  std::array<std::complex<double>, Dim * Dim> matR{};
  for (std::size_t i = 0; i < Dim; i++) {
    for (std::size_t j = 0; j < Dim; j++) {
      for (std::size_t k = 0; k < Dim; k++) {
        matR[i + Dim * j] += m[i + Dim * k] * eigenpair.vectors[k + Dim * j];
      }
    }
  }
  // compute R*L
  std::array<std::complex<double>, Dim * Dim> rL{};
  for (std::size_t i = 0; i < Dim; i++) {
    for (std::size_t j = 0; j < Dim; j++) {
      rL[i + Dim * j] = eigenpair.vectors[i + Dim * j] * eigenpair.values[j];
    }
  }
  // compare residual
  for (std::size_t i = 0; i < Dim * Dim; i++) {
    REQUIRE(std::abs(matR[i] - rL[i]) == AbsApprox(0.0).epsilon(Epsilon));
  }
}

TEST_CASE("Eigenvalues are correctly computed" * doctest::test_suite("numerical")) {
  SUBCASE("Eigenvalues 3") {
    constexpr std::size_t Dim = 3;

    std::array<std::array<std::complex<double>, Dim * Dim>, 2> matrices = {{
        {2.0, -3.0, -3.0, -2.0, 1.0, -2.0, -1.0, 1.0, 4.0},
        {-2.0, 3.0, 3.0, 2.0, -1.0, 2.0, 1.0, -1.0, -4.0},
    }};
    for (auto& m : matrices) {
      seissol::eigenvalues::Eigenpair<std::complex<double>, Dim> eigenpair{};
      seissol::eigenvalues::computeEigenvalues(m, eigenpair);
      testResidual<Dim>(m, eigenpair);
    }
  }
}

TEST_CASE("The eigenvectors have unit length and a real, positive largest component" *
          doctest::test_suite("numerical")) {
  constexpr std::size_t Dim = 3;
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 1e2;

  // a real matrix with the simple real eigenvalues -1, 3 and 5, whose eigenvectors are then real,
  // and a complex one
  std::array<std::array<std::complex<double>, Dim * Dim>, 2> matrices = {{
      {2.0, -3.0, -3.0, -2.0, 1.0, -2.0, -1.0, 1.0, 4.0},
      {std::complex<double>(1.0, 2.0),
       std::complex<double>(0.5, -1.0),
       std::complex<double>(0.0, 3.0),
       std::complex<double>(2.0, 0.0),
       std::complex<double>(-1.0, 1.0),
       std::complex<double>(0.25, 0.5),
       std::complex<double>(1.0, -2.0),
       std::complex<double>(3.0, 1.0),
       std::complex<double>(-2.0, 0.0)},
  }};
  for (std::size_t index = 0; index < matrices.size(); ++index) {
    seissol::eigenvalues::Eigenpair<std::complex<double>, Dim> eigenpair{};
    seissol::eigenvalues::computeEigenvalues(matrices[index], eigenpair);
    testResidual<Dim>(matrices[index], eigenpair);
    for (std::size_t j = 0; j < Dim; ++j) {
      double norm = 0.0;
      std::size_t largest = 0;
      for (std::size_t i = 0; i < Dim; ++i) {
        const auto& component = eigenpair.vectors[i + Dim * j];
        norm += std::norm(component);
        if (std::abs(component) > std::abs(eigenpair.vectors[largest + Dim * j])) {
          largest = i;
        }
        if (index == 0) {
          REQUIRE(component.imag() == AbsApprox(0.0).epsilon(Epsilon));
        }
      }
      REQUIRE(std::sqrt(norm) == AbsApprox(1.0).epsilon(Epsilon));
      REQUIRE(eigenpair.vectors[largest + Dim * j].imag() == 0.0);
      REQUIRE(eigenpair.vectors[largest + Dim * j].real() > 0.0);
    }
  }
}

/// Checks the eigenpairs of k A, the operator of a planar wave along x with the wavenumber k in a
/// solid of the density rho and the stiffness given as the upper triangle of its Voigt matrix
/// (xx, yy, zz, yz, xz, xy), for the quantities (sxx, syy, szz, sxy, syz, sxz, u, v, w). Its rows
/// of the stresses hold stiffnesses, its rows of the velocities the inverse of the density. Its
/// eigenvalues are -k c, 0 three times and k c, with the speeds c of the Christoffel matrix of x,
/// and the velocities of each wave lie in the eigenspace of that matrix that its speed belongs to.
void testPlaneWaveOperator(const std::array<double, 21>& upper, double rho, double k) {
  constexpr std::size_t Dim = 9;
  std::array<std::array<double, 6>, 6> stiffness{};
  std::size_t next = 0;
  for (std::size_t i = 0; i < 6; ++i) {
    for (std::size_t j = i; j < 6; ++j) {
      stiffness[i][j] = upper[next];
      stiffness[j][i] = upper[next];
      ++next;
    }
  }
  // the rows of the stresses in the Voigt matrix, and the ones of the strain rates that the
  // velocities u, v and w drive along x
  constexpr std::array<std::size_t, 6> StressRow = {0, 1, 2, 5, 3, 4};
  constexpr std::array<std::size_t, 3> StrainRow = {0, 5, 4};

  std::array<std::complex<double>, Dim * Dim> m{};
  for (std::size_t s = 0; s < 6; ++s) {
    for (std::size_t v = 0; v < 3; ++v) {
      m[s + Dim * (6 + v)] = -k * stiffness[StressRow[s]][StrainRow[v]];
    }
  }
  m[6 + Dim * 0] = -k / rho;
  m[7 + Dim * 3] = -k / rho;
  m[8 + Dim * 5] = -k / rho;

  Eigen::Matrix3d christoffel;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      christoffel(i, j) = stiffness[StrainRow[i]][StrainRow[j]] / rho;
    }
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> waves(christoffel);
  const double largest = k * std::sqrt(waves.eigenvalues()(2));
  constexpr double Tolerance = 1e-12;

  seissol::eigenvalues::Eigenpair<std::complex<double>, Dim> eigenpair{};
  seissol::eigenvalues::computeEigenvalues(m, eigenpair);

  // sorted by their real parts: the waves to -x from the fastest, the null space, the waves to +x
  // from the slowest
  for (std::size_t wave = 0; wave < 3; ++wave) {
    const double frequency = k * std::sqrt(waves.eigenvalues()(wave));
    // the projection onto the polarizations of the waves of this speed
    Eigen::Matrix3cd polarization = Eigen::Matrix3cd::Zero();
    for (std::size_t other = 0; other < 3; ++other) {
      if (std::abs(k * std::sqrt(waves.eigenvalues()(other)) - frequency) <= Tolerance * largest) {
        const Eigen::Vector3cd direction =
            waves.eigenvectors().col(other).cast<std::complex<double>>();
        polarization += direction * direction.adjoint();
      }
    }
    for (const auto& [index, sign] : {std::pair{2 - wave, -1.0}, std::pair{6 + wave, 1.0}}) {
      REQUIRE(eigenpair.values[index].real() ==
              AbsApprox(sign * frequency).epsilon(Tolerance * largest));
      REQUIRE(eigenpair.values[index].imag() == AbsApprox(0.0).epsilon(Tolerance * largest));

      Eigen::Vector3cd velocity;
      for (std::size_t i = 0; i < 3; ++i) {
        velocity(i) = eigenpair.vectors[6 + i + Dim * index];
      }
      const double across = (velocity - polarization * velocity).norm() / velocity.norm();
      REQUIRE(across == AbsApprox(0.0).epsilon(Tolerance));
    }
  }
  for (std::size_t index = 3; index < 6; ++index) {
    REQUIRE(std::abs(eigenpair.values[index]) == AbsApprox(0.0).epsilon(Tolerance * largest));
  }
}

TEST_CASE("Eigenvalues of a badly scaled operator are correctly computed" *
          doctest::test_suite("numerical")) {
  // a transversely isotropic solid with axes tilted 40 degrees from z, in which all three waves
  // along x couple all components of the velocity
  const std::array<double, 21> tilted = {
      9.267e10, 2.631e10,  2.722e10,  1.407e09, -5.096e09, -3.703e09, // c11 to c16
      8.515e10, 3.032e10,  -6.695e09, 2.043e09, -2.809e09,            // c22 to c26
      7.83e10,  -3.837e09, -2.216e09, 2.684e09,                       // c33 to c36
      2.964e10, 3.938e08,  -1.108e09,                                 // c44 to c46
      2.919e10, -4.051e09,                                            // c55 and c56
      3.182e10};                                                      // c66
  // an isotropic solid with a P wave of 6000 m/s and two S waves of 3464 m/s
  const std::array<double, 21> isotropic = {9.72e10, 3.24e10, 3.24e10, 0.0, 0.0, 0.0, // c11 to c16
                                            9.72e10, 3.24e10, 0.0,     0.0, 0.0,      // c22 to c26
                                            9.72e10, 0.0,     0.0,     0.0,           // c33 to c36
                                            3.24e10, 0.0,     0.0,                    // c44 to c46
                                            3.24e10, 0.0,                             // c55 and c56
                                            3.24e10};                                 // c66
  constexpr double Rho = 2700.0;
  constexpr double K = 2 * std::numbers::pi / 100;
  for (const auto& stiffness : {tilted, isotropic}) {
    testPlaneWaveOperator(stiffness, Rho, K);
  }
}

TEST_CASE("Balancing leaves entries that are not finite alone" * doctest::test_suite("numerical")) {
  // a density of zero puts an infinite entry into the rows of the velocities
  Eigen::Matrix<std::complex<double>, 3, 3> m;
  m << 1.0, std::numeric_limits<double>::infinity(), 0.0, 1e10, 2.0, 1e-4, 0.0, 1e-4, 3.0;
  const auto scale = seissol::eigenvalues::balance(m);
  for (Eigen::Index i = 0; i < scale.size(); ++i) {
    REQUIRE(std::isfinite(scale(i)));
  }
}

} // namespace seissol::unit_test
