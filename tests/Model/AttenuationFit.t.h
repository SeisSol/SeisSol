// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Equations/viscoacoustic/Model/Attenuation.h"
#include "Equations/viscoelastic/Model/Attenuation.h"

#include <complex>
#include <cstddef>

namespace seissol::unit_test {

namespace {

constexpr double FreqCentral = 1.0;
constexpr double FreqRatio = 100.0;
/// Frequencies inside the fitted band. A fit with a handful of mechanisms
/// ripples around the target across the band, so the tolerance is generous;
/// it is still tight enough to catch a wrong sign, a wrong modulus, or a
/// factor.
constexpr double Frequencies[] = {0.1, 0.3, 1.0, 3.0, 10.0};
constexpr double Tolerance = 0.1;

/// The quality factor the fitted parameters actually produce at `frequency`,
/// from the standard relaxation sum: M(w) = M_u - sum_l dM_l w_l / (w_l + i w).
template <typename MaterialT, typename Unrelaxed, typename Defect>
double fittedQ(const MaterialT& material, double frequency, Unrelaxed unrelaxed, Defect defect) {
  const double angular = 2 * M_PI * frequency;
  std::complex<double> modulus(unrelaxed(material), 0.0);
  for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
    modulus -= defect(material, mech) * material.omega[mech] /
               std::complex<double>(material.omega[mech], angular);
  }
  return std::abs(modulus.real() / modulus.imag());
}

} // namespace

TEST_CASE("Attenuation fit reproduces the requested Q" * doctest::test_suite("model")) {
  SUBCASE("viscoelastic, both wave types") {
    model::ViscoElasticMaterial<3> material;
    material.rho = 2700.0;
    material.mu = 2.0e10;
    material.lambda = 3.0e10;
    material.qp = 60.0;
    material.qs = 30.0;
    const double requestedP = material.qp;
    const double requestedS = material.qs;

    physics::fitAttenuation(material, FreqCentral, FreqRatio);

    for (const double frequency : Frequencies) {
      const double qp = fittedQ(
          material,
          frequency,
          [](const auto& m) { return m.lambda + 2 * m.mu; },
          [](const auto& m, std::size_t l) { return m.getDeltaLambda(l) + 2 * m.getDeltaMu(l); });
      const double qs = fittedQ(
          material,
          frequency,
          [](const auto& m) { return m.mu; },
          [](const auto& m, std::size_t l) { return m.getDeltaMu(l); });
      CHECK(qp == doctest::Approx(requestedP).epsilon(Tolerance));
      CHECK(qs == doctest::Approx(requestedS).epsilon(Tolerance));
    }
  }

  SUBCASE("viscoacoustic, the one wave type there is") {
    model::ViscoAcousticMaterial<3> material;
    material.rho = 1000.0;
    material.lambda = 2.2e9;
    material.qp = 40.0;
    const double requested = material.qp;

    physics::fitAttenuation(material, FreqCentral, FreqRatio);

    for (const double frequency : Frequencies) {
      // lambda is the bulk modulus here, so the bulk defect is the whole
      // defect; dividing it by three, as the accessor once did, shows up as a
      // Q three times too large.
      const double q = fittedQ(
          material,
          frequency,
          [](const auto& m) { return m.lambda; },
          [](const auto& m, std::size_t l) { return m.getDeltaBulk(l); });
      CHECK(q == doctest::Approx(requested).epsilon(Tolerance));
    }
  }
}

TEST_CASE("Attenuation fit leaves a well posed material" * doctest::test_suite("model")) {
  SUBCASE("viscoelastic") {
    model::ViscoElasticMaterial<3> material;
    material.rho = 2700.0;
    material.mu = 2.0e10;
    material.lambda = 3.0e10;
    material.qp = 60.0;
    material.qs = 30.0;
    physics::fitAttenuation(material, FreqCentral, FreqRatio);

    // The relaxed moduli are what a static load sees, so they have to be
    // positive and below the unrelaxed ones.
    CHECK(material.getMuRelaxed() > 0.0);
    CHECK(material.getMuRelaxed() < material.mu);
    CHECK(material.getBulkRelaxed() > 0.0);
    CHECK(material.attenuationWellPosed());
  }

  SUBCASE("viscoacoustic") {
    model::ViscoAcousticMaterial<3> material;
    material.rho = 1000.0;
    material.lambda = 2.2e9;
    material.qp = 40.0;
    physics::fitAttenuation(material, FreqCentral, FreqRatio);

    CHECK(material.getBulkRelaxed() > 0.0);
    CHECK(material.getBulkRelaxed() < material.lambda);
    CHECK(material.attenuationWellPosed());
  }
}

} // namespace seissol::unit_test
