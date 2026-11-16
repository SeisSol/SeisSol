// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/InitProcedure/Internal/MeshLayout.h"
#include "Initializer/InitProcedure/Internal/NeighborLinks.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Tree/Colormap.h"

#include <cstddef>
#include <map>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::initializer;

namespace meshlayouttest {

inline LTSColorMap layoutColorMap() {
  return LTSColorMap(EnumLayer<HaloType>({HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
                     EnumLayer<std::size_t>({0, 1}),
                     EnumLayer<ConfigId>({0}));
}

} // namespace meshlayouttest

TEST_CASE("Combining the copy and interior cells of a cluster" *
          doctest::test_suite("initializer")) {
  const auto colorMap = meshlayouttest::layoutColorMap();
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
    const auto combined =
        seissol::initializer::internal::combineCopyInterior(layout, colorMap, {6});
    CHECK(combined.moved == 1);
    CHECK(combined.copyLayers[color(HaloType::Copy, 0)]);
    CHECK_FALSE(combined.copyLayers[color(HaloType::Copy, 1)]);

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
    CHECK(seissol::initializer::internal::combineCopyInterior(layout, colorMap, {0}).moved == 0);

    CHECK(layout[color(HaloType::Copy, 0)].cellMap.size() == 2);
    CHECK(layout[color(HaloType::Interior, 0)].cellMap.size() == 3);
  }
}

TEST_CASE("Combining the copy and interior cells with a threshold per configuration" *
          doctest::test_suite("initializer")) {
  // e.g. orders 6 and 3 in one run: a layer of the cheaper configuration costs less per cell, while
  // its fixed part is much the same, so it combines up to more cells
  const LTSColorMap colorMap(
      EnumLayer<HaloType>({HaloType::Ghost, HaloType::Copy, HaloType::Interior}),
      EnumLayer<std::size_t>({0}),
      EnumLayer<ConfigId>({0, 1}));
  const auto color = [&](HaloType halo, ConfigId config) {
    return colorMap.color(halo, 0, config);
  };

  std::vector<ClusterMap> layout(colorMap.size());
  for (const ConfigId config : {ConfigId{0}, ConfigId{1}}) {
    layout[color(HaloType::Copy, config)].cellMap = {0, 1};
    layout[color(HaloType::Interior, config)].cellMap = {2, 3, 4};
  }

  const auto combined =
      seissol::initializer::internal::combineCopyInterior(layout, colorMap, {4, 6});

  CHECK(combined.moved == 1);
  CHECK(layout[color(HaloType::Interior, 0)].cellMap.size() == 3);
  CHECK(layout[color(HaloType::Copy, 1)].cellMap.size() == 5);
  CHECK(layout[color(HaloType::Interior, 1)].cellMap.empty());
}

TEST_CASE("Combining the copy and interior fault faces of a cluster" *
          doctest::test_suite("initializer")) {
  const auto colorMap = meshlayouttest::layoutColorMap();
  const auto color = [&](HaloType halo, std::size_t lts) {
    return colorMap.color(halo, lts, ConfigId{0});
  };

  std::vector<std::vector<std::size_t>> layout(colorMap.size());
  layout[color(HaloType::Copy, 0)] = {3};
  layout[color(HaloType::Interior, 0)] = {1, 2};
  layout[color(HaloType::Copy, 1)] = {4, 5, 6};
  layout[color(HaloType::Interior, 1)] = {7, 8, 9};

  const auto below = seissol::initializer::internal::faceLayersBelow(layout, colorMap, 4);
  CHECK(below[color(HaloType::Copy, 0)]);
  CHECK_FALSE(below[color(HaloType::Copy, 1)]);
  CHECK_FALSE(below[color(HaloType::Interior, 0)]);
  CHECK(seissol::initializer::internal::combineCopyInteriorDR(layout, colorMap, below) == 1);

  // the faces of the copy layer go to the interior layer, which then also waits for the ghost cells
  CHECK(layout[color(HaloType::Copy, 0)].empty());
  CHECK(layout[color(HaloType::Interior, 0)] == std::vector<std::size_t>{1, 2, 3});
  CHECK(layout[color(HaloType::Copy, 1)].size() == 3);
  CHECK(layout[color(HaloType::Interior, 1)].size() == 3);
}

TEST_CASE("Combining the copy and interior cells where a choice per cluster says so" *
          doctest::test_suite("initializer")) {
  const auto colorMap = meshlayouttest::layoutColorMap();
  const auto color = [&](HaloType halo, std::size_t lts) {
    return colorMap.color(halo, lts, ConfigId{0});
  };

  std::vector<ClusterMap> layout(colorMap.size());
  layout[color(HaloType::Copy, 0)].cellMap = {0, 1};
  layout[color(HaloType::Copy, 0)].regions.emplace_back(7, 0, 0, 2, 1);
  layout[color(HaloType::Interior, 0)].cellMap = {2, 3, 4};
  layout[color(HaloType::Copy, 1)].cellMap = {5, 6, 7, 8};
  layout[color(HaloType::Copy, 1)].regions.emplace_back(8, 0, 0, 4, 1);
  layout[color(HaloType::Interior, 1)].cellMap = {9, 10, 11, 12, 13, 14};

  // the larger cluster only, which the choice sees with the sizes of its two layers
  const auto combined = seissol::initializer::internal::combineCopyInteriorWhere(
      layout,
      colorMap,
      [&](std::size_t copyColor, std::size_t copyCells, std::size_t interiorCells) {
        return copyColor == color(HaloType::Copy, 1) && copyCells == 4 && interiorCells == 6;
      });
  CHECK(combined.moved == 1);
  CHECK_FALSE(combined.copyLayers[color(HaloType::Copy, 0)]);
  CHECK(combined.copyLayers[color(HaloType::Copy, 1)]);
  CHECK(layout[color(HaloType::Interior, 0)].cellMap.size() == 3);
  CHECK(layout[color(HaloType::Copy, 1)].cellMap.size() == 10);
  CHECK(layout[color(HaloType::Interior, 1)].cellMap.empty());
}

TEST_CASE("The time the halo of a copy layer takes" * doctest::test_suite("initializer")) {
  const auto colorMap = meshlayouttest::layoutColorMap();
  const auto color = [&](HaloType halo, std::size_t lts) {
    return colorMap.color(halo, lts, ConfigId{0});
  };

  std::vector<ClusterMap> layout(colorMap.size());
  // the copy cells of cluster 0 wait for 10 ghost cells from rank 1 and 5 from rank 2; those of
  // cluster 1 for 20 from rank 1
  layout[color(HaloType::Ghost, 0)].regions.emplace_back(1, 0, color(HaloType::Copy, 0), 10, 1);
  layout[color(HaloType::Ghost, 0)].regions.emplace_back(2, 0, color(HaloType::Copy, 0), 5, 2);
  layout[color(HaloType::Ghost, 1)].regions.emplace_back(3, 0, color(HaloType::Copy, 1), 20, 1);

  std::map<int, seissol::initializer::internal::NeighborLink> links;
  links[1] = {.latency = 10e-6, .bytesPerSecond = 1e9};
  links[2] = {.latency = 50e-6};
  const std::vector<std::size_t> bytesPerCell(colorMap.size(), 1000);

  const auto halo =
      seissol::initializer::internal::haloTimes(layout, colorMap, links, bytesPerCell);
  // in microseconds (doctest::Approx compares small numbers with a margin of about 1e-5 otherwise)
  // the slowest of the ranks: rank 2 by its latency
  CHECK(halo[color(HaloType::Copy, 0)] * 1e6 == doctest::Approx(50));
  // 20000 bytes at 1 GB/s after 10 us
  CHECK(halo[color(HaloType::Copy, 1)] * 1e6 == doctest::Approx(30));
  CHECK(halo[color(HaloType::Interior, 0)] == 0.0);
}

} // namespace seissol::unit_test
