// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/viscoelastic/Model/Datastructures.h"
#include "Expr/SderivFrontend.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Model/MaterialType.h"
#include "Physics/InstantaneousTimeMirrorManager.h"
#include "Reader/Scripting/CompiledReader.h"
#include "TestHelper.h"

#include <array>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Anisotropic Instantaneous Time Mirror supports only BothWaves reflection type" *
          doctest::skip(!std::is_same_v<model::MaterialT, model::AnisotropicMaterial>) *
          doctest::test_suite("physics")) {
  using seissol::initializer::parameters::ReflectionType;

  CHECK(seissol::physics::isAnisotropicReflectionTypeSupported(ReflectionType::BothWaves));
  CHECK_FALSE(seissol::physics::isAnisotropicReflectionTypeSupported(ReflectionType::Pwave));
  CHECK_FALSE(seissol::physics::isAnisotropicReflectionTypeSupported(ReflectionType::Swave));
  CHECK_FALSE(
      seissol::physics::isAnisotropicReflectionTypeSupported(ReflectionType::BothWavesVelocity));
}

TEST_CASE("Instantaneous Time Mirror Swave lambda scaling matches implementation formula" *
          doctest::test_suite("physics")) {
  const double lambda = 8.0;
  const double mu = 3.0;
  const double velocityScalingFactor = 2.0;
  const double expected =
      (lambda + 2.0 * mu) / velocityScalingFactor - 2.0 * velocityScalingFactor * mu;

  CHECK(seissol::physics::getSwaveScaledLambda(lambda, mu, velocityScalingFactor) ==
        AbsApprox(expected).epsilon(1.0e-15));
}

TEST_CASE("Elastic Instantaneous Time Mirror time-step scaling depends on reflection type" *
          doctest::test_suite("physics")) {
  using seissol::initializer::parameters::ReflectionType;

  const double velocityScalingFactor = 2.5;
  CHECK(seissol::physics::getElasticTimeStepScalingFactor(ReflectionType::BothWaves,
                                                          velocityScalingFactor) ==
        AbsApprox(1.0 / velocityScalingFactor).epsilon(1.0e-15));
  CHECK(seissol::physics::getElasticTimeStepScalingFactor(ReflectionType::Pwave,
                                                          velocityScalingFactor) ==
        AbsApprox(1.0 / velocityScalingFactor).epsilon(1.0e-15));
  CHECK(seissol::physics::getElasticTimeStepScalingFactor(ReflectionType::Swave,
                                                          velocityScalingFactor) ==
        AbsApprox(velocityScalingFactor).epsilon(1.0e-15));
  CHECK(seissol::physics::getElasticTimeStepScalingFactor(ReflectionType::BothWavesVelocity,
                                                          velocityScalingFactor) ==
        AbsApprox(1.0).epsilon(1.0e-15));
}
namespace itm {

/// The materials `materials` after the mirror script `source` with the factor `n`.
template <typename MaterialT>
void mirror(const std::string& source,
            double n,
            std::vector<MaterialT>& materials,
            const initializer::parameters::ModelParameters& parameters = {}) {
  reader::scripting::CompiledReader reader(expr::compileSderivModule(source), nullptr);
  physics::checkMirrorScript<MaterialT>(reader, "test");
  std::vector<MaterialT*> pointers;
  std::vector<std::array<double, 3>> centers;
  for (auto& material : materials) {
    pointers.push_back(&material);
    centers.push_back({0.5, -0.25, 1.0});
  }
  physics::mirrorMaterials<MaterialT>(reader, "test", n, pointers, centers, parameters);
}

} // namespace itm

TEST_CASE("Instantaneous Time Mirror: the scripts of the reflection types change an elastic "
          "material as the built-in mirror did" *
          doctest::test_suite("physics")) {
  using seissol::initializer::parameters::ReflectionType;
  // small enough that the S wave mirror keeps lambda positive
  const double n = 1.2;
  const double rho = 2.7;
  const double mu = 3.1;
  const double lambda = 8.3;
  for (const auto type : {ReflectionType::BothWaves,
                          ReflectionType::BothWavesVelocity,
                          ReflectionType::Pwave,
                          ReflectionType::Swave}) {
    CAPTURE(static_cast<int>(type));
    const auto source = physics::mirrorScript(model::MaterialType::Elastic, type);
    REQUIRE(source.has_value());
    std::vector<model::ElasticMaterial> materials(1, model::ElasticMaterial({rho, mu, lambda}));
    itm::mirror(*source, n, materials);
    const auto& material = materials.front();
    // the scaling of the built-in mirror, in the same order of operations
    switch (type) {
    case ReflectionType::BothWaves:
      CHECK(material.rho == rho);
      CHECK(material.mu == mu * n * n);
      CHECK(material.lambda == lambda * n * n);
      break;
    case ReflectionType::BothWavesVelocity:
      CHECK(material.rho == rho * n);
      CHECK(material.mu == mu * n);
      CHECK(material.lambda == lambda * n);
      break;
    case ReflectionType::Pwave:
      CHECK(material.rho == rho);
      CHECK(material.mu == mu);
      CHECK(material.lambda == lambda * n * n);
      break;
    case ReflectionType::Swave:
      CHECK(material.rho == n * rho);
      CHECK(material.mu == n * mu);
      CHECK(material.lambda ==
            AbsApprox(physics::getSwaveScaledLambda(lambda, mu, n)).epsilon(1.0e-15));
      break;
    }
  }
}

TEST_CASE("Instantaneous Time Mirror: the mirror of an anisotropic material scales its density" *
          doctest::test_suite("physics")) {
  using seissol::initializer::parameters::ReflectionType;
  CHECK_FALSE(physics::mirrorScript(model::MaterialType::Anisotropic, ReflectionType::Pwave));
  CHECK_FALSE(physics::mirrorScript(model::MaterialType::Acoustic, ReflectionType::BothWaves));
  const auto source =
      physics::mirrorScript(model::MaterialType::Anisotropic, ReflectionType::BothWaves);
  REQUIRE(source.has_value());
  const model::ElasticMaterial elastic({2.0, 3.0, 4.0});
  std::vector<model::AnisotropicMaterial> materials(1, model::AnisotropicMaterial(elastic));
  const double c11 = materials.front().c11;
  itm::mirror(*source, 1.5, materials);
  CHECK(materials.front().rho == 2.0 / (1.5 * 1.5));
  CHECK(materials.front().c11 == c11);
}

TEST_CASE("Instantaneous Time Mirror: a script changes any parameter of a material" *
          doctest::test_suite("physics")) {
  SUBCASE("viscoelastic: Q and the moduli of the material file, fitted anew") {
    initializer::parameters::ModelParameters parameters;
    parameters.freqCentral = 0.3;
    parameters.freqRatio = 100;
    model::ViscoElasticMaterial<3> material;
    material.rho = 2787.0;
    material.mu = 3.0e10;
    material.lambda = 3.1e10;
    material.qp = 58.0;
    material.qs = 29.0;
    material.initialize(parameters);
    std::vector<model::ViscoElasticMaterial<3>> materials(1, material);
    itm::mirror("out def Qp = Qp0 * n\nout def mu = mu0 * n * n\n", 2.0, materials, parameters);

    // the same as a material with these parameters from the start
    model::ViscoElasticMaterial<3> expected;
    expected.rho = 2787.0;
    expected.mu = 3.0e10 * 2.0 * 2.0;
    expected.lambda = 3.1e10;
    expected.qp = 116.0;
    expected.qs = 29.0;
    expected.initialize(parameters);
    const auto& mirrored = materials.front();
    CHECK(mirrored.qp == 116.0);
    CHECK(mirrored.qs == 29.0);
    CHECK(mirrored.mu == expected.mu);
    CHECK(mirrored.lambda == expected.lambda);
    std::array<double, 5> given{};
    mirrored.scriptParameters(given.data());
    CHECK(given[1] == 1.2e11);
    CHECK(given[2] == 3.1e10);
  }

  SUBCASE("poroelastic") {
    model::PoroElasticMaterial material;
    material.rho = 2500.0;
    material.porosity = 0.2;
    material.permeability = 1e-12;
    std::vector<model::PoroElasticMaterial> materials(1, material);
    itm::mirror("out def porosity = porosity0 * 0.5\nout def rho = rho0 + x\n", 2.0, materials);
    CHECK(materials.front().porosity == 0.1);
    CHECK(materials.front().rho == 2500.5);
    CHECK(materials.front().permeability == 1e-12);
  }

  SUBCASE("the inputs of a mirror") {
    CHECK(physics::mirrorInputs({"rho", "Qp"}) ==
          std::vector<std::string>{"x", "y", "z", "n", "rho0", "Qp0"});
  }
}

} // namespace seissol::unit_test
