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

TEST_CASE("Backmap") {
  constexpr std::size_t DupCount = 4;

  std::vector<std::size_t> dataToMesh{0, 1, 2, 3, 4, 3, 6, 4, 9, 8, 6, 2, 3, 5, 0, 2};
  std::vector<std::size_t> layers{5, 6, 5};

  std::unordered_map<std::size_t, std::vector<initializer::StoragePosition>> compare;

  auto backmap = initializer::StorageBackmap<DupCount>();

  const std::size_t cellCount = 12;

  backmap.setSize(cellCount);

  // set up backmap and comparison

  std::size_t global = 0;
  const auto* zero = dataToMesh.data();
  for (const auto [layerId, layerSize] : common::enumerate(layers)) {
    const auto* zeroLayer = dataToMesh.data() + global;
    for (std::size_t i = 0; i < layerSize; ++i) {
      backmap.addElement(layerId, zero, zeroLayer, dataToMesh[global], i);

      compare[dataToMesh[global]].emplace_back(layerId, i, global);

      ++global;
    }
  }

  // compare

  for (std::size_t i = 0; i < cellCount; ++i) {
    for (std::size_t j = 0; j < DupCount; ++j) {
      const auto pos = backmap.getDup(i, j);
      if (j < compare[i].size()) {
        REQUIRE(pos.has_value());
        REQUIRE(pos.value() == compare[i][j]);
      } else {
        REQUIRE(!pos.has_value());
      }
    }

    if (!compare[i].empty()) {
      REQUIRE(backmap.get(i) == compare[i][0]);
    }
  }
}

} // namespace seissol::unit_test
