// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigRegistry.h"
#include "Memory/Tree/Colormap.h"

#include <cstddef>
#include <vector>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Colormap" * doctest::test_suite("memory")) {
  const initializer::LTSColorMap colorMap(initializer::EnumLayer(std::vector<HaloType>{
                                              HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
                                          initializer::EnumLayer(std::vector<std::size_t>{1, 2, 3}),
                                          initializer::EnumLayer(std::vector<ConfigId>{0, 1}));

  REQUIRE(colorMap.size() == 18);
  CHECK(colorMap.argument(0).lts == 1);
  CHECK(colorMap.argument(0).halo == HaloType::Interior);
  CHECK(colorMap.argument(0).config == 0);
  CHECK(colorMap.argument(1).lts == 1);
  CHECK(colorMap.argument(1).halo == HaloType::Copy);
  CHECK(colorMap.argument(1).config == 0);
  CHECK(colorMap.argument(9).lts == 1);
  CHECK(colorMap.argument(9).halo == HaloType::Interior);
  CHECK(colorMap.argument(9).config == 1);

  for (std::size_t color = 0; color < colorMap.size(); ++color) {
    CHECK(colorMap.colorId(colorMap.argument(color)) == color);
  }
}

} // namespace seissol::unit_test
