// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Model/Quantities.h"

#include <cstddef>
#include <vector>

namespace seissol::unit_test {

namespace {
using ViscoAcoustic = model::ViscoAcousticMaterial<3>;

/// A material with the moduli already fitted, i.e. lambda unrelaxed and theta
/// holding the negated defects.
ViscoAcoustic fitted() {
  ViscoAcoustic material;
  material.rho = 1000.0;
  material.lambda = 2.4e9;
  for (std::size_t mech = 0; mech < ViscoAcoustic::Mechanisms; ++mech) {
    material.omega[mech] = 0.6 * std::pow(10.0, static_cast<double>(mech));
    material.theta[mech][0] = -1.0e8 * (mech + 1);
  }
  return material;
}
} // namespace

TEST_CASE("Visco-acoustic moduli" * doctest::test_suite("equations")) {
  const auto material = fitted();

  SUBCASE("the bulk defect is the whole defect") {
    // An acoustic material has no shear modulus, so lambda is the bulk
    // modulus. The elastic formula, dLambda + 2/3 dMu, does not carry over:
    // taking a third of the defect makes the branch energies and the
    // dissipation rate three times too small.
    for (std::size_t mech = 0; mech < ViscoAcoustic::Mechanisms; ++mech) {
      CHECK(material.getDeltaBulk(mech) == material.getDeltaLambda(mech));
      CHECK(material.getDeltaLambda(mech) == -material.theta[mech][0]);
      CHECK(material.getDeltaBulk(mech) > 0.0);
    }
  }

  SUBCASE("relaxed is unrelaxed minus every defect") {
    double expected = material.lambda;
    for (std::size_t mech = 0; mech < ViscoAcoustic::Mechanisms; ++mech) {
      expected -= material.getDeltaBulk(mech);
    }
    CHECK(material.getBulkRelaxed() == doctest::Approx(expected));
    CHECK(material.getBulkRelaxed() == material.getLambdaRelaxed());
    CHECK(material.getLambdaUnrelaxed() == material.lambda);
    CHECK(material.getBulkRelaxed() < material.getLambdaUnrelaxed());
  }

  SUBCASE("a material that relaxes below zero is rejected") {
    auto broken = material;
    broken.theta[0][0] = -2 * broken.lambda;
    CHECK_FALSE(broken.attenuationWellPosed());
    CHECK(material.attenuationWellPosed());
  }
}

TEST_CASE("Visco-acoustic quantity layout" * doctest::test_suite("equations")) {
  SUBCASE("pressure is a scalar and the memory variables are too") {
    // The elastic layout has a symmetric tensor here. Carrying that over is
    // what once made the pressure rotate with two velocity components.
    constexpr auto Traction =
        model::roleKind(ViscoAcoustic::PrimaryGroups, model::FaceRole::Traction);
    static_assert(Traction == model::QuantityKind::Scalar);
    static_assert(ViscoAcoustic::MechanismGroups.size() == 1);
    static_assert(ViscoAcoustic::MechanismGroups[0].kind == model::QuantityKind::Scalar);
    static_assert(ViscoAcoustic::NumberPerMechanism == 1);
  }

  SUBCASE("the velocities start after the single pressure") {
    static_assert(ViscoAcoustic::VelocityOffset == 1);
    static_assert(ViscoAcoustic::TractionComponents == 1);
    static_assert(ViscoAcoustic::NumElasticQuantities == 4);
    static_assert(ViscoAcoustic::NumQuantities == 4 + ViscoAcoustic::Mechanisms);
  }

  SUBCASE("one memory variable per mechanism, each on its own") {
    // Three mechanisms are three one-component blocks, not one three-component
    // block: consecutive mechanisms overlapped when that was confused.
    constexpr auto Rotation = ViscoAcoustic::RotationGroups;
    std::size_t scalars = 0;
    for (const auto& group : Rotation) {
      scalars += static_cast<std::size_t>(group.kind == model::QuantityKind::Scalar);
    }
    // the pressure plus one per mechanism, in the layout that carries them
    const std::size_t expected =
        1 + (ViscoAcoustic::RotationRepetitions * ViscoAcoustic::MechanismGroups.size());
    CHECK(scalars == expected);
  }
}

TEST_CASE("Visco-acoustic vector constructor" * doctest::test_suite("equations")) {
  // Two values per mechanism, a relaxation frequency and the one source entry
  // an acoustic mechanism has. Reading three of them, as the elastic layout
  // would, runs into the next mechanism.
  std::vector<double> values{1000.0, 2.4e9};
  for (std::size_t mech = 0; mech < ViscoAcoustic::Mechanisms; ++mech) {
    values.push_back(10.0 + mech);
    values.push_back(-1.0e8 * (mech + 1));
  }
  REQUIRE(values.size() == ViscoAcoustic::Parameters);

  const ViscoAcoustic material(values);
  CHECK(material.rho == 1000.0);
  CHECK(material.lambda == 2.4e9);
  for (std::size_t mech = 0; mech < ViscoAcoustic::Mechanisms; ++mech) {
    CHECK(material.omega[mech] == 10.0 + mech);
    CHECK(material.theta[mech][0] == -1.0e8 * (mech + 1));
  }
}

} // namespace seissol::unit_test
