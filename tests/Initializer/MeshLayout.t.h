// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/InitProcedure/Internal/MeshLayout.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Tree/Colormap.h"

#include <cstddef>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::initializer;

namespace {

LTSColorMap layoutColorMap() {
  return LTSColorMap(EnumLayer<HaloType>({HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
                     EnumLayer<std::size_t>({0, 1}),
                     EnumLayer<ConfigId>({0}));
}

} // namespace

TEST_CASE("Combining the copy and interior cells of a cluster" *
          doctest::test_suite("initializer")) {
  const auto colorMap = layoutColorMap();
  const auto color = [&](HaloType halo, std::size_t lts) {
    return colorMap.color(halo, lts, ConfigId{0});
  };

  std::vector<ClusterMap> layout(colorMap.size());
  // cluster 0: 2 copy cells sent to rank 1, 3 interior cells; cluster 1: 4 copy cells, 6 interior
  layout[color(HaloType::Copy, 0)].cellMap = {0, 1};
  layout[color(HaloType::Copy, 0)].regions.emplace_back(7, 0, 0, 2, 1);
  layout[color(HaloType::Interior, 0)].cellMap = {2, 3, 4};
  layout[color(HaloType::Copy, 1)].cellMap = {5, 6, 7, 8};
  layout[color(HaloType::Copy, 1)].regions.emplace_back(8, 0, 0, 4, 1);
  layout[color(HaloType::Interior, 1)].cellMap = {9, 10, 11, 12, 13, 14};

  SUBCASE("below the threshold, the interior cells follow the copy cells in the copy layer") {
    CHECK(seissol::initializer::internal::combineCopyInterior(layout, colorMap, 6) == 1);

    CHECK(layout[color(HaloType::Copy, 0)].cellMap == std::vector<std::size_t>{0, 1, 2, 3, 4});
    CHECK(layout[color(HaloType::Interior, 0)].cellMap.empty());
    // the halo stays as it was: the regions cover the leading copy cells
    REQUIRE(layout[color(HaloType::Copy, 0)].regions.size() == 1);
    CHECK(layout[color(HaloType::Copy, 0)].regions[0].count == 2);

    // 10 cells are not below 6
    CHECK(layout[color(HaloType::Copy, 1)].cellMap.size() == 4);
    CHECK(layout[color(HaloType::Interior, 1)].cellMap.size() == 6);
  }

  SUBCASE("a threshold of 0 combines nothing") {
    CHECK(seissol::initializer::internal::combineCopyInterior(layout, colorMap, 0) == 0);

    CHECK(layout[color(HaloType::Copy, 0)].cellMap.size() == 2);
    CHECK(layout[color(HaloType::Interior, 0)].cellMap.size() == 3);
  }
}

TEST_CASE("Combining the copy and interior fault faces of a cluster" *
          doctest::test_suite("initializer")) {
  const auto colorMap = layoutColorMap();
  const auto color = [&](HaloType halo, std::size_t lts) {
    return colorMap.color(halo, lts, ConfigId{0});
  };

  std::vector<std::vector<std::size_t>> layout(colorMap.size());
  layout[color(HaloType::Copy, 0)] = {3};
  layout[color(HaloType::Interior, 0)] = {1, 2};
  layout[color(HaloType::Copy, 1)] = {4, 5, 6};
  layout[color(HaloType::Interior, 1)] = {7, 8, 9};

  CHECK(seissol::initializer::internal::combineCopyInteriorDR(layout, colorMap, 4) == 1);

  // the faces of the copy layer go to the interior layer, which then also waits for the ghost cells
  CHECK(layout[color(HaloType::Copy, 0)].empty());
  CHECK(layout[color(HaloType::Interior, 0)] == std::vector<std::size_t>{1, 2, 3});
  CHECK(layout[color(HaloType::Copy, 1)].size() == 3);
  CHECK(layout[color(HaloType::Interior, 1)].size() == 3);
}

} // namespace seissol::unit_test
