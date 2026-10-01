// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Real.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Solver/MultipleSimulations.h"

#include <cstddef>
#include <string>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Built configurations" * doctest::test_suite("common")) {
  REQUIRE(builtConfigCount() > 0);

  const auto id = findConfig(Config::Value);

  SUBCASE("The configuration of the generated code is built") {
    REQUIRE(id.has_value());
    CHECK(configValue(id.value()) == Config::Value);
    CHECK(findConfig(configName(Config::Value)) == id);
    CHECK(defaultConfig() < builtConfigCount());
  }

  SUBCASE("Nothing else is found") {
    auto other = Config::Value;
    other.convergenceOrder += 10;
    CHECK_FALSE(findConfig(other).has_value());
    CHECK_FALSE(findConfig(configName(other)).has_value());
    CHECK_FALSE(findConfig(std::string("not-a-configuration")).has_value());
  }

  SUBCASE("The layout matches the generated code") {
    REQUIRE(id.has_value());
    const auto& layout = configLayout(id.value());
    const auto order = Config::ConvergenceOrder;
    CHECK(layout.config == Config::Value);
    CHECK(layout.numQuantities == model::MaterialT::NumQuantities);
    CHECK(layout.numBasisFunctions == order * (order + 1) * (order + 2) / 6);
    CHECK(layout.basisFunctionDimension == multisim::BasisFunctionDimension);
    CHECK(layout.dofsSize == tensor::Q::size());
    // (the unknowns of a cell may hold fewer quantities than the material, e.g. when a solver
    // stores the anelastic ones apart)
    const auto storedQuantities = tensor::Q::Shape[multisim::BasisFunctionDimension + 1];
    CHECK(storedQuantities <= layout.numQuantities);
    CHECK(layout.dofsSize >= layout.numBasisFunctions * storedQuantities * Config::NumSimulations);
    CHECK(layout.drNumPoints == dr::misc::NumBoundaryGaussPoints);
    CHECK(layout.drNumPaddedPoints == dr::misc::NumPaddedPoints);
    CHECK(layout.drNumPaddedPoints >= layout.drNumPoints * Config::NumSimulations);
    CHECK(layout.drNumQuantities == dr::misc::NumQuantities);
    CHECK(layout.drNumTimePoints == dr::misc::TimeSteps);
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

} // namespace seissol::unit_test
