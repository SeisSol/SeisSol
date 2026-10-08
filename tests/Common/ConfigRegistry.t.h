// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Solver/MultipleSimulations.h"
#include "TestHelper.h"

#include <cctype>
#include <cstddef>
#include <optional>
#include <string>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Built configurations" * doctest::test_suite("common")) {
  REQUIRE(builtConfigCount() > 0);
  CHECK(defaultConfig() < builtConfigCount());

  SUBCASE("Every configuration of the generated code is built under its id") {
    forEachConfig([](auto cfg) {
      using Cfg = decltype(cfg);
      const auto id = configIdOf<Cfg>();
      CAPTURE(id);
      CHECK(configValue(id) == Cfg::Value);
      CHECK(findConfig(Cfg::Value) == id);
      CHECK(findConfig(configName(Cfg::Value)) == id);
    });
  }

  SUBCASE("Nothing else is found") {
    auto other = configValue(defaultConfig());
    other.convergenceOrder += 10;
    CHECK_FALSE(findConfig(other).has_value());
    CHECK_FALSE(findConfig(configName(other)).has_value());
    CHECK_FALSE(findConfig(std::string("not-a-configuration")).has_value());
  }

  SUBCASE("The layout matches the generated code") {
    forEachConfig([](auto cfg) {
      using Cfg = decltype(cfg);
      CAPTURE(configIdOf<Cfg>());
      const auto& layout = configLayout(configIdOf<Cfg>());
      const auto order = Cfg::ConvergenceOrder;
      CHECK(layout.config == Cfg::Value);
      CHECK(layout.numQuantities == model::MaterialOf<Cfg>::NumQuantities);
      CHECK(layout.velocityOffset == model::MaterialOf<Cfg>::VelocityOffset);
      CHECK(layout.numBasisFunctions == order * (order + 1) * (order + 2) / 6);
      CHECK(layout.basisFunctionDimension == multisim::BasisDim<Cfg>);
      CHECK(layout.dofsSize == tensor::Q<Cfg>::size());
      // (the unknowns of a cell may hold fewer quantities than the material, e.g. when a solver
      // stores the anelastic ones apart)
      const auto storedQuantities = tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
      CHECK(storedQuantities <= layout.numQuantities);
      CHECK(layout.dofsSize >= layout.numBasisFunctions * storedQuantities * Cfg::NumSimulations);
      CHECK(layout.drNumPoints == dr::misc::NumBoundaryGaussPoints<Cfg>);
      CHECK(layout.drNumPaddedPoints == dr::misc::NumPaddedPoints<Cfg>);
      CHECK(layout.drNumPaddedPoints >= layout.drNumPoints * Cfg::NumSimulations);
      CHECK(layout.drNumQuantities == dr::misc::NumQuantities<Cfg>);
      CHECK(layout.drNumTimePoints == dr::misc::TimeSteps<Cfg>);
    });
  }

  SUBCASE("Every layout agrees with the kernels of its id in runtime.h") {
    for (std::size_t other = 0; other < builtConfigCount(); ++other) {
      const auto config = static_cast<ConfigId>(other);
      CAPTURE(config);
      const auto& layout = configLayout(config);

      const auto* dofs = runtime::init::Q::descriptor(config);
      REQUIRE(dofs != nullptr);
      CHECK(yateto::sizeOf(dofs->datatype) == sizeOfRealType(layout.config.precision));
      CHECK(dofs->size == layout.dofsSize);
      CHECK(dofs->shape[layout.basisFunctionDimension] == layout.numBasisFunctions);

      const auto* face = runtime::init::QInterpolated::descriptor(config);
      REQUIRE(face != nullptr);
      CHECK(face->shape[layout.basisFunctionDimension] == layout.drNumPoints);
      CHECK(face->shape[layout.basisFunctionDimension + 1] == layout.drNumQuantities);
      CHECK(face->size == layout.drNumPaddedPoints * layout.drNumQuantities);
    }
  }

  SUBCASE("Every built configuration is described") {
    const auto description = describeBuiltConfigs();
    for (std::size_t other = 0; other < builtConfigCount(); ++other) {
      CHECK(description.find(configName(configValue(static_cast<ConfigId>(other)))) !=
            std::string::npos);
    }
  }
}

TEST_CASE(
    "The default configuration is the one SEISSOL_CONFIGURATION names and otherwise the first" *
    doctest::test_suite("common")) {
  const auto last = static_cast<ConfigId>(builtConfigCount() - 1);

  SUBCASE("named, spelled as in a parameter file") {
    auto name = configName(configValue(last));
    for (auto& character : name) {
      character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
    }
    const ScopedEnvironment environment("SEISSOL_CONFIGURATION", " " + name + " ");
    CHECK(environmentConfig() == last);
    CHECK(defaultConfig() == last);
  }

  SUBCASE("not set") {
    const ScopedEnvironment environment("SEISSOL_CONFIGURATION", std::nullopt);
    CHECK_FALSE(environmentConfig().has_value());
    CHECK(defaultConfig() == 0);
  }
}

} // namespace seissol::unit_test
