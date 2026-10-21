// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Iterator.h"
#include "Memory/Tree/Backmap.h"
namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Backmap" * doctest::test_suite("memory")) {
  constexpr std::size_t DupCount = 4;

  std::vector<std::size_t> dataToMesh{0, 1, 2, 3, 4, 3, 6, 4, 9, 8, 6, 2, 3, 5, 0, 2};
  std::vector<std::size_t> layers{5, 6, 5};

  std::unordered_map<std::size_t, std::vector<initializer::StoragePosition>> compare;

  auto backmap = initializer::StorageBackmap<DupCount>();

  const std::size_t cellCount = 12;

  backmap.setSize(cellCount, 1);

  // set up backmap and comparison

  std::size_t global = 0;
  const auto* zero = dataToMesh.data();
  for (const auto [layerId, layerSize] : common::enumerate(layers)) {
    const auto* zeroLayer = dataToMesh.data() + global;
    for (std::size_t i = 0; i < layerSize; ++i) {
      backmap.addElement(layerId, zero, zeroLayer, dataToMesh[global], 0, i);

      compare[dataToMesh[global]].emplace_back(layerId, i, global);

      ++global;
    }
  }

  // compare

  for (std::size_t i = 0; i < cellCount; ++i) {
    for (std::size_t j = 0; j < DupCount; ++j) {
      const auto pos = backmap.getDup(i, 0, j);
      if (j < compare[i].size()) {
        REQUIRE(pos.has_value());
        CHECK(pos.value() == compare[i][j]);
      } else {
        REQUIRE(!pos.has_value());
      }
    }

    if (!compare[i].empty()) {
      CHECK(backmap.get(i, 0) == compare[i][0]);
    }
  }
}

TEST_CASE("Backmap with copies" * doctest::test_suite("memory")) {
  constexpr std::size_t DupCount = 2;
  constexpr std::size_t CellCount = 5;
  constexpr std::size_t CopyCount = 3;

  auto backmap = initializer::StorageBackmap<DupCount>();
  backmap.setSize(CellCount, CopyCount);
  CHECK(backmap.copyCount() == CopyCount);

  // a first layer with the copies of every cell next to each other, as the layout has them; and a
  // second one with a duplicate of one copy of a cell
  std::vector<std::size_t> data(CellCount * CopyCount + 1);
  const auto* zero = data.data();
  for (std::size_t cell = 0; cell < CellCount; ++cell) {
    for (std::size_t copy = 0; copy < CopyCount; ++copy) {
      CHECK(backmap.addElement(0, zero, zero, cell, copy, cell * CopyCount + copy) == 0);
    }
  }
  const auto* secondLayer = zero + CellCount * CopyCount;
  CHECK(backmap.addElement(1, zero, secondLayer, 2, 1, 0) == 1);

  for (std::size_t cell = 0; cell < CellCount; ++cell) {
    for (std::size_t copy = 0; copy < CopyCount; ++copy) {
      const auto index = cell * CopyCount + copy;
      const auto position = backmap.get(cell, copy);
      CHECK(position == initializer::StoragePosition(0, index, index));
      CHECK(position.global == index);

      const auto duplicate = backmap.getDup(cell, copy, 1);
      if (cell == 2 && copy == 1) {
        REQUIRE(duplicate.has_value());
        CHECK(duplicate.value() == initializer::StoragePosition(1, 0, CellCount * CopyCount));
        CHECK(duplicate.value().global == CellCount * CopyCount);
      } else {
        // a duplicate belongs to its copy only
        CHECK(!duplicate.has_value());
      }
    }
  }
}

} // namespace seissol::unit_test
