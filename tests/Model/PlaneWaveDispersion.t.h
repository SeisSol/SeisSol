// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_PLANEWAVEDISPERSION_T_H_
#define SEISSOL_TESTS_MODEL_PLANEWAVEDISPERSION_T_H_

#include <doctest.h>

#include "Core/Common/MaterialType.h"
#include "Core/Numerical/Eigenvalues.h"
#include "TestConfigs.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Physics/viscoelastic/Attenuation.h"
#include "WavePropagation/Solver/MaterialOperators.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

// The plane-wave operator k A - i E of a material, whose eigenvalues are the frequencies omega of
// its plane waves exp(i (omega t - k x)) (see physics::Planarwave), against what the theory of the
// material says they are. The planar waves of the convergence cases take their exact solution from
// the same operator, so that only a comparison with the theory holds the coefficient matrices and
// the source terms themselves.

namespace seissol::unit_test {

namespace planewaves {

using Complex = std::complex<double>;

/// The eigenvalues of the plane-wave operator of `material` along x, for the wave number `k`.
template <typename Cfg>
std::vector<Complex> frequencies(const model::MaterialOf<Cfg>& material, double k) {
  constexpr auto Count = model::MaterialOf<Cfg>::NumQuantities;
  const double direction[3] = {k, 0.0, 0.0};
  std::array<Complex, Count * Count> planeWaveOperator{};
  model::getPlaneWaveOperator<Cfg>(material, direction, planeWaveOperator.data());
  eigenvalues::Eigenpair<Complex, Count> eigenpair;
  eigenvalues::computeEigenvalues(planeWaveOperator, eigenpair);
  return {eigenpair.values.begin(), eigenpair.values.end()};
}

/// The roots of the polynomial sum_i coefficients[i] x^i, from the eigenvalues of its companion
/// matrix in the variable x / scale, whose coefficients are of a size alike.
inline std::vector<Complex> roots(const std::vector<Complex>& coefficients, double scale) {
  const auto degree = static_cast<Eigen::Index>(coefficients.size()) - 1;
  std::vector<Complex> scaled(coefficients.size());
  for (std::size_t i = 0; i < coefficients.size(); ++i) {
    scaled[i] = coefficients[i] * std::pow(scale, static_cast<double>(i));
  }
  Eigen::MatrixXcd companion = Eigen::MatrixXcd::Zero(degree, degree);
  for (Eigen::Index i = 0; i < degree; ++i) {
    companion(0, i) = -scaled[degree - 1 - i] / scaled[degree];
    if (i + 1 < degree) {
      companion(i + 1, i) = 1.0;
    }
  }
  const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(companion);
  std::vector<Complex> found;
  for (const auto& root : solver.eigenvalues()) {
    found.push_back(root * scale);
  }
  return found;
}

/// The distance of `value` to the nearest of `values`.
inline double distanceTo(const std::vector<Complex>& values, Complex value) {
  double nearest = std::numeric_limits<double>::infinity();
  for (const auto& candidate : values) {
    nearest = std::min(nearest, std::abs(candidate - value));
  }
  return nearest;
}

} // namespace planewaves

TEST_CASE_TEMPLATE_DEFINE("The plane waves of a viscoelastic material run at the speed and decay "
                          "with the quality factor it is given" *
                              doctest::test_suite("model"),
                          Cfg,
                          ViscoelasticPlaneWaves) {
  using namespace planewaves;

  // The moduli are those of the wave speeds at the central frequency of the band (see
  // fitAttenuation), and the quality factors those the mechanisms fit over the band. With an
  // exp(i (omega t - k x)) of a real k, rho omega^2 = k^2 M(omega) gives a frequency of the real
  // part k c and the quality factor Re omega / (2 Im omega), up to terms of 1 / Q^2.
  constexpr double FreqCentral = 1.0;
  constexpr double FreqRatio = 100.0;
  constexpr double Rho = 2700.0;
  constexpr double Mu = 3.24e10;
  constexpr double Lambda = 3.24e10;
  constexpr double Qp = 40.0;
  constexpr double Qs = 20.0;

  model::MaterialOf<Cfg> material;
  material.rho = Rho;
  material.mu = Mu;
  material.lambda = Lambda;
  material.qp = Qp;
  material.qs = Qs;
  physics::fitAttenuation(material, FreqCentral, FreqRatio);

  const double central = 2 * std::numbers::pi * FreqCentral;
  // (not as a structured binding: CAPTURE takes the speed into a lambda, which clang 18 does not
  // allow for a structured binding under OpenMP)
  for (const auto& mode :
       {std::pair{std::sqrt((Lambda + 2 * Mu) / Rho), Qp}, std::pair{std::sqrt(Mu / Rho), Qs}}) {
    const double speed = mode.first;
    const double quality = mode.second;
    CAPTURE(speed);
    const auto omegas = frequencies<Cfg>(material, central / speed);

    // no mode grows
    for (const auto& omega : omegas) {
      CHECK(omega.imag() >= -1e-9 * central);
    }

    // the wave, of the frequency next to the central one
    const auto wave = *std::min_element(omegas.begin(), omegas.end(), [&](auto a, auto b) {
      return std::abs(a - central) < std::abs(b - central);
    });
    CHECK(wave.real() == doctest::Approx(central).epsilon(1 / (2 * quality * quality)));
    CHECK(wave.real() / (2 * wave.imag()) == doctest::Approx(quality).epsilon(0.1));

    // The planar wave takes its modes by their places in this order (physics::Planarwave): P
    // last, and second an S wave that travels the other way, at the frequency -conj(omega) of
    // the one that travels along x. The modes that decay without travelling sort between.
    if (quality == Qp) {
      CHECK(omegas.back() == wave);
    } else {
      CHECK(std::abs(omegas[1] + std::conj(wave)) <= 1e-9 * central);
    }
  }
}

TEST_CASE_TEMPLATE_DEFINE("The plane waves of a poroelastic material solve Biot's dispersion "
                          "relation" *
                              doctest::test_suite("model"),
                          Cfg,
                          PoroelasticPlaneWaves) {
  using namespace planewaves;

  // A sandstone saturated with water. In the notation of docs/physical-models.rst and of Carcione
  // (Wave Fields in Real Media, 7.2 and 7.4), with u the displacement of the solid and w the one of
  // the fluid relative to it, times the porosity:
  //   rho u'' + rho_F w'' = div sigma,  sigma = lambda_M tr(eps) I + 2 mu_M eps - alpha p I,
  //   rho_F u'' + m w'' + b w' = -grad p,  p = -M (alpha div u + div w),
  // with the density rho = (1 - phi) rho_S + phi rho_F, the mass coupling m = T rho_F / phi, the
  // friction b = nu / kappa, alpha = 1 - K_M / K_S for the bulk modulus K_M of the matrix, and
  // M = K_S / (1 - phi - K_M / K_S + phi K_S / K_F).
  model::MaterialOf<Cfg> material;
  material.bulkSolid = 40.0e9;
  material.rho = 2500.0;
  material.lambda = 12.0e9;
  material.mu = 10.0e9;
  material.porosity = 0.2;
  material.permeability = 6.0e-13;
  material.tortuosity = 3.0;
  material.bulkFluid = 2.5e9;
  material.rhoFluid = 1040.0;
  material.viscosity = 0.001;

  const double phi = material.porosity;
  const double rho = (1 - phi) * material.rho + phi * material.rhoFluid;
  const double rhoF = material.rhoFluid;
  const double m = material.tortuosity * material.rhoFluid / phi;
  const double b = material.viscosity / material.permeability;
  const double bulkMatrix = material.lambda + 2 * material.mu / 3;
  const double alpha = 1 - bulkMatrix / material.bulkSolid;
  const double modulus = material.bulkSolid / (1 - phi - bulkMatrix / material.bulkSolid +
                                               phi * material.bulkSolid / material.bulkFluid);
  const double h = material.lambda + 2 * material.mu + alpha * alpha * modulus;
  const double mu = material.mu;

  // waves of a kilometre and of a centimetre: the slow wave diffuses at the first, and propagates
  // at the second, above Biot's frequency b phi / (2 pi T rho_F), some 50 kHz
  for (const double wavelength : {1000.0, 0.01}) {
    CAPTURE(wavelength);
    const double k = 2 * std::numbers::pi / wavelength;
    const double k2 = k * k;
    const Complex i(0.0, 1.0);

    // exp(i (omega t - k x)) along x. The P waves: det [[k^2 H - omega^2 rho, k^2 alpha M -
    // omega^2 rho_F], [k^2 alpha M - omega^2 rho_F, k^2 M - omega^2 m + i omega b]] = 0
    const double a = k2 * h;
    const double c = k2 * modulus;
    const double d = k2 * alpha * modulus;
    const double scale = k * std::sqrt(h / rho);
    const std::vector<Complex> pWaves = roots({Complex(a * c - d * d),
                                               i * a * b,
                                               Complex(-a * m - rho * c + 2 * d * rhoF),
                                               -i * rho * b,
                                               Complex(rho * m - rhoF * rhoF)},
                                              scale);
    // each S wave: (k^2 mu - omega^2 rho)(i b - omega m) - omega^3 rho_F^2 = 0, one for each
    // direction of polarization
    const std::vector<Complex> sWaves = roots(
        {i * b * k2 * mu, Complex(-k2 * mu * m), -i * rho * b, Complex(rho * m - rhoF * rhoF)},
        scale);

    const auto omegas = frequencies<Cfg>(material, k);
    REQUIRE(omegas.size() == 13);

    std::vector<Complex> expected(pWaves);
    expected.insert(expected.end(), sWaves.begin(), sWaves.end());
    expected.insert(expected.end(), sWaves.begin(), sWaves.end());
    // the stresses that a wave along x leaves alone
    expected.insert(expected.end(), 3, Complex(0.0, 0.0));
    REQUIRE(expected.size() == omegas.size());

    // the slow wave diffuses at frequencies far above those of the others
    const auto tolerance = [&](Complex omega) { return 1e-6 * std::max(scale, std::abs(omega)); };
    for (const auto& omega : expected) {
      CAPTURE(omega);
      CHECK(distanceTo(omegas, omega) <= tolerance(omega));
    }
    for (const auto& omega : omegas) {
      CAPTURE(omega);
      CHECK(distanceTo(expected, omega) <= tolerance(omega));
      // no mode grows
      CHECK(omega.imag() >= -tolerance(omega));
    }
  }
}

TEST_CASE_TEMPLATE_APPLY(ViscoelasticPlaneWaves,
                         ConfigsOfMaterial<model::MaterialType::Viscoelastic>);
TEST_CASE_TEMPLATE_APPLY(PoroelasticPlaneWaves,
                         ConfigsOfMaterial<model::MaterialType::Poroelastic>);

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_PLANEWAVEDISPERSION_T_H_
