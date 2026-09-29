// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/BasicTypedefs.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"

#include <cstddef>
#include <vector>

namespace seissol::unit_test {

using namespace seissol;

struct TestDescriptor {
  struct Var1 : public initializer::Variable<int> {};
  struct Var2 : public initializer::Variable<float> {};
  struct Var3 : public initializer::Variable<double[123]> {};
  struct Bucket : public initializer::Bucket<float> {};
  struct Scratchpad : public initializer::Scratchpad<float> {};
};

TEST_CASE("Storage" * doctest::test_suite("memory")) {
  initializer::Storage<initializer::GenericVarmap> storage;

  // NOTE: the LTSColorMap is hard-coded to the storage right now.
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2, 3}),
      initializer::TraitLayer(std::vector<initializer::ConfigVariant>{Config()}));

  constexpr auto Alignment1 = sizeof(void*);
  constexpr auto Alignment2 = sizeof(void*) * 4;

  storage.add<TestDescriptor::Var1>(Ghost, Alignment1, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var2>(Copy, Alignment2, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var3>(Interior, Alignment2, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Bucket>(
      initializer::LayerMask(), Alignment1, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Scratchpad>(
      initializer::LayerMask(), Alignment1, initializer::AllocationMode::HostOnly);

  storage.setLayerCount(colorMap);

  REQUIRE(storage.numChildren() == colorMap.size());

  storage.fixate();

  for (const auto [i, layer] : common::enumerate(storage.leaves())) {
    CHECK(layer.getIdentifier().lts == colorMap.argument(i).lts);
    CHECK(layer.getIdentifier().halo == colorMap.argument(i).halo);
    CHECK(layer.getIdentifier().config.index() == colorMap.argument(i).config.index());
  }

  for (auto [i, layer] : common::enumerate(storage.leaves())) {
    layer.setNumberOfCells(i + 1);
  }

  storage.allocateVariables();
  storage.touchVariables();

  for (const auto [i, layer] : common::enumerate(storage.leaves())) {
    CHECK(layer.size() == i + 1);
  }

  storage.allocateBuckets();
  storage.allocateScratchPads();
}

TEST_CASE("Storage scratchpad sharing" * doctest::test_suite("memory")) {
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2}),
      initializer::TraitLayer(std::vector<initializer::ConfigVariant>{Config()}));

  const auto setup = [&](initializer::Storage<initializer::GenericVarmap>& storage) {
    storage.add<TestDescriptor::Scratchpad>(
        initializer::LayerMask(), sizeof(void*), initializer::AllocationMode::HostOnly);
    storage.setLayerCount(colorMap);
    storage.fixate();
    for (auto [i, layer] : common::enumerate(storage.leaves())) {
      layer.setNumberOfCells(1);
      // one layer without any demand
      layer.setEntrySize<TestDescriptor::Scratchpad>(i == 2 ? 0 : (i + 1) * sizeof(float) * 16);
    }
    storage.allocateVariables();
    storage.allocateBuckets();
  };

  SUBCASE("Shared by all layers") {
    initializer::Storage<initializer::GenericVarmap> storage;
    setup(storage);
    storage.allocateScratchPads(initializer::ScratchpadSharing::Shared);

    const float* first = nullptr;
    for (auto& layer : storage.leaves()) {
      const auto* scratchpad = layer.var<TestDescriptor::Scratchpad>();
      REQUIRE(scratchpad != nullptr);
      if (first == nullptr) {
        first = scratchpad;
      }
      CHECK(scratchpad == first);
    }
  }

  SUBCASE("One per layer") {
    initializer::Storage<initializer::GenericVarmap> storage;
    setup(storage);
    storage.allocateScratchPads(initializer::ScratchpadSharing::PerLayer);

    std::vector<const float*> scratchpads;
    for (auto [i, layer] : common::enumerate(storage.leaves())) {
      const auto* scratchpad = layer.var<TestDescriptor::Scratchpad>();
      if (i == 2) {
        CHECK(scratchpad == nullptr);
      } else {
        REQUIRE(scratchpad != nullptr);
        scratchpads.push_back(scratchpad);
      }
    }
    for (std::size_t i = 0; i < scratchpads.size(); ++i) {
      for (std::size_t j = i + 1; j < scratchpads.size(); ++j) {
        CHECK(scratchpads[i] != scratchpads[j]);
      }
    }

    // the demands must not overlap either
    for (auto [i, layer] : common::enumerate(storage.leaves())) {
      auto* scratchpad = layer.var<TestDescriptor::Scratchpad>();
      if (scratchpad != nullptr) {
        const auto count = layer.getEntrySize<TestDescriptor::Scratchpad>() / sizeof(float);
        for (std::size_t k = 0; k < count; ++k) {
          scratchpad[k] = static_cast<float>(i);
        }
      }
    }
    for (auto [i, layer] : common::enumerate(storage.leaves())) {
      const auto* scratchpad = layer.var<TestDescriptor::Scratchpad>();
      if (scratchpad != nullptr) {
        const auto count = layer.getEntrySize<TestDescriptor::Scratchpad>() / sizeof(float);
        for (std::size_t k = 0; k < count; ++k) {
          CHECK(scratchpad[k] == static_cast<float>(i));
        }
      }
    }
  }
}

} // namespace seissol::unit_test
