// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/acoustic/Model/Datastructures.h"
#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/viscoacoustic/Model/Datastructures.h"
#include "Equations/viscoelastic/Model/Datastructures.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Model/CommonDatastructures.h"
#include "Model/ScriptParameters.h"

#include <cstddef>
#include <string>
#include <vector>

namespace seissol::unit_test::script_parameters {

/// The parameters of `material` as scripts read them, through the base class.
template <typename MaterialT>
std::vector<double> parametersOf(const MaterialT& material) {
  std::vector<double> values(model::ScriptParameterCount<MaterialT>, -1.0);
  static_cast<const model::Material&>(material).scriptParameters(values.data());
  return values;
}

TEST_CASE("ScriptParameters: the parameters of a material in order, as the material file gave "
          "them" *
          doctest::test_suite("model")) {
  using model::ScriptParameterCount;
  using model::scriptParameterNames;
  using Names = std::vector<std::string>;

  SUBCASE("elastic") {
    CHECK(scriptParameterNames<model::ElasticMaterial>() == Names{"rho", "mu", "lambda"});
    CHECK(ScriptParameterCount<model::ElasticMaterial> == 3);
    model::ElasticMaterial material;
    material.rho = 2.0;
    material.mu = 3.0;
    material.lambda = 4.0;
    CHECK(parametersOf(material) == std::vector<double>{2.0, 3.0, 4.0});
  }

  SUBCASE("acoustic: a fluid has no shear modulus") {
    CHECK(scriptParameterNames<model::AcousticMaterial>() == Names{"rho", "lambda", "mu"});
    model::AcousticMaterial material;
    material.rho = 1.5;
    material.lambda = 2.25;
    CHECK(parametersOf(material) == std::vector<double>{1.5, 2.25, 0.0});
  }

  SUBCASE("anisotropic: mu and lambda of the closest isotropic material") {
    const auto names = scriptParameterNames<model::AnisotropicMaterial>();
    REQUIRE(names.size() == 24);
    CHECK(ScriptParameterCount<model::AnisotropicMaterial> == 24);
    CHECK(names[0] == "rho");
    CHECK(names[1] == "c11");
    CHECK(names[21] == "c66");
    CHECK(names[22] == "mu");
    CHECK(names[23] == "lambda");
    model::ElasticMaterial elastic;
    elastic.rho = 2.0;
    elastic.mu = 3.0;
    elastic.lambda = 4.0;
    const model::AnisotropicMaterial material(elastic);
    const auto values = parametersOf(material);
    CHECK(values[0] == 2.0);
    // c11 = lambda + 2 mu, c12 = lambda, c44 = mu
    CHECK(values[1] == 10.0);
    CHECK(values[2] == 4.0);
    CHECK(values[16] == 3.0);
    CHECK(values[22] == doctest::Approx(3.0));
    CHECK(values[23] == doctest::Approx(4.0));
  }

  SUBCASE("poroelastic") {
    CHECK(scriptParameterNames<model::PoroElasticMaterial>() == Names{"rho",
                                                                      "mu",
                                                                      "lambda",
                                                                      "bulk_solid",
                                                                      "porosity",
                                                                      "permeability",
                                                                      "tortuosity",
                                                                      "bulk_fluid",
                                                                      "rho_fluid",
                                                                      "viscosity"});
    model::PoroElasticMaterial material;
    material.rho = 1.0;
    material.mu = 2.0;
    material.lambda = 3.0;
    material.bulkSolid = 4.0;
    material.porosity = 0.5;
    material.permeability = 6.0;
    material.tortuosity = 7.0;
    material.bulkFluid = 8.0;
    material.rhoFluid = 9.0;
    material.viscosity = 10.0;
    CHECK(parametersOf(material) ==
          std::vector<double>{1.0, 2.0, 3.0, 4.0, 0.5, 6.0, 7.0, 8.0, 9.0, 10.0});
  }

  SUBCASE("viscoelastic: the moduli the material file gave, not those of the fit") {
    CHECK(scriptParameterNames<model::ViscoElasticMaterial<3>>() ==
          Names{"rho", "mu", "lambda", "Qp", "Qs"});
    model::ViscoElasticMaterial<3> material;
    material.rho = 2787.0;
    material.mu = 3.0e10;
    material.lambda = 3.1e10;
    material.qp = 58.0;
    material.qs = 29.0;
    initializer::parameters::ModelParameters parameters;
    parameters.freqCentral = 0.3;
    parameters.freqRatio = 100;
    material.initialize(parameters);
    // the fit replaces the moduli by the unrelaxed ones
    CHECK(material.mu != 3.0e10);
    CHECK(parametersOf(material) == std::vector<double>{2787.0, 3.0e10, 3.1e10, 58.0, 29.0});
  }

  SUBCASE("viscoacoustic") {
    CHECK(scriptParameterNames<model::ViscoAcousticMaterial<2>>() ==
          Names{"rho", "lambda", "Qp", "mu"});
    model::ViscoAcousticMaterial<2> material;
    material.rho = 1000.0;
    material.lambda = 2.25e9;
    material.qp = 40.0;
    initializer::parameters::ModelParameters parameters;
    parameters.freqCentral = 1.0;
    parameters.freqRatio = 100;
    material.initialize(parameters);
    CHECK(material.lambda != 2.25e9);
    CHECK(parametersOf(material) == std::vector<double>{1000.0, 2.25e9, 40.0, 0.0});
  }

  SUBCASE("a name, by its position") {
    CHECK(model::scriptParameterIndex<model::ElasticMaterial>("lambda") == 2);
    CHECK(!model::scriptParameterIndex<model::ElasticMaterial>("porosity").has_value());
  }
}

} // namespace seissol::unit_test::script_parameters
