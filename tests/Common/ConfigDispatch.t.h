// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"

#include <cstddef>
#include <variant>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Configuration dispatch" * doctest::test_suite("common")) {
  REQUIRE(builtConfigCount() == std::variant_size_v<ConfigVariant>);

  SUBCASE("Every id leads to the type of its configuration") {
    for (std::size_t id = 0; id < builtConfigCount(); ++id) {
      const auto value = dispatchConfig(static_cast<ConfigId>(id),
                                        [](auto config) { return decltype(config)::Value; });
      CHECK(value == configValue(static_cast<ConfigId>(id)));
    }
  }

  SUBCASE("Every type leads back to its id") {
    forEachConfig([](auto config) {
      using Cfg = decltype(config);
      static_assert(configIdOf<Cfg>() < std::variant_size_v<ConfigVariant>);
      CHECK(findConfig(Cfg::Value) == configIdOf<Cfg>());
    });
  }

  SUBCASE("An id outside the build is rejected") {
    CHECK_THROWS(dispatchConfig(static_cast<ConfigId>(builtConfigCount()),
                                [](auto /*config*/) { return 0; }));
  }
}

} // namespace seissol::unit_test
