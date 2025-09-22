// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <Common/ConfigHelper.h>
#include <Memory/Tree/Colormap.h>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Colormap") {
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2, 3}),
      initializer::TraitLayer(std::vector<ConfigVariant>{ConfigVariantList[0]}));

  REQUIRE(colorMap.size() == 9);
  REQUIRE(colorMap.argument(0).lts == 1);
  REQUIRE(colorMap.argument(0).halo == HaloType::Interior);
  REQUIRE(colorMap.argument(0).config.index() == 0);
  REQUIRE(colorMap.argument(1).lts == 1);
  REQUIRE(colorMap.argument(1).halo == HaloType::Copy);
  REQUIRE(colorMap.argument(1).config.index() == 0);
}

} // namespace seissol::unit_test
