// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Config.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Iterator.h"
#include "Core/Common/Real.h"
#include "Core/Common/TypeId.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Tree/Colormap.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "Initializer/BasicTypedefs.h"
#include "TestConfigs.h"

#include <cstddef>
#include <type_traits>
#include <typeinfo>
#include <vector>

namespace seissol::unit_test {

using namespace seissol;

template <typename Cfg>
using PerConfigArray = RealT<Cfg::Precision>[Cfg::ConvergenceOrder + 1];

struct TestDescriptor {
  struct Var1 : public initializer::Variable<int> {};
  struct Var2 : public initializer::Variable<float> {};
  struct Var3 : public initializer::Variable<double[123]> {};
  struct Var4 : public initializer::VariantVariable<PerConfigArray> {};
  struct Bucket : public initializer::Bucket<float> {};
  struct Scratchpad : public initializer::Scratchpad<float> {};
};

TEST_CASE_TEMPLATE("Storage" * doctest::test_suite("memory"), Cfg, SEISSOL_CONFIG_TYPES) {
  initializer::Storage storage;

  // NOTE: the LTSColorMap is hard-coded to the storage right now.
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2, 3}),
      initializer::EnumLayer(std::vector<ConfigId>{configIdOf<Cfg>()}));

  constexpr auto Alignment1 = sizeof(void*);
  constexpr auto Alignment2 = sizeof(void*) * 4;

  storage.add<TestDescriptor::Var1>(Ghost, Alignment1, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var2>(Copy, Alignment2, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var3>(Interior, Alignment2, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var4>(
      initializer::LayerMask(), Alignment1, initializer::AllocationMode::HostOnly);
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
    CHECK(layer.getIdentifier().config == colorMap.argument(i).config);
  }

  // a variable that depends on the configuration is sized by the configuration of its layer
  for (const auto& layer : storage.leaves()) {
    CHECK(storage.info<TestDescriptor::Var4>().bytesLayer(layer.getIdentifier()) ==
          sizeof(PerConfigArray<Cfg>));
  }

  for (auto [i, layer] : common::enumerate(storage.leaves())) {
    layer.setNumberOfCells(i + 1);
  }

  storage.allocateVariables();
  storage.touchVariables();

  for (const auto [i, layer] : common::enumerate(storage.leaves())) {
    CHECK(layer.size() == i + 1);
  }

  for (auto& layer : storage.leaves()) {
    auto* perConfig = layer.var<TestDescriptor::Var4>(Cfg());
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t j = 0; j <= Cfg::ConvergenceOrder; ++j) {
        CHECK(perConfig[cell][j] == 0);
        perConfig[cell][j] = static_cast<RealT<Cfg::Precision>>(cell + j);
      }
    }
  }

  // a cell of a configuration sees the variables as that configuration holds them, whichever way
  // it is reached
  for (auto [color, layer] : common::enumerate(storage.leaves())) {
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      auto ref = layer.cellRef<Cfg>(cell);
      static_assert(
          std::is_same_v<decltype(ref.template get<TestDescriptor::Var4>()), PerConfigArray<Cfg>&>);
      static_assert(
          std::is_same_v<decltype(ref.template get<TestDescriptor::Var3>()), double (&)[123]>);
      CHECK(&ref.template get<TestDescriptor::Var4>() ==
            &layer.var<TestDescriptor::Var4>(Cfg())[cell]);

      initializer::StoragePosition position;
      position.color = color;
      position.cell = cell;
      CHECK(&storage.lookup<TestDescriptor::Var4>(Cfg(), position) ==
            &ref.template get<TestDescriptor::Var4>());
      CHECK(&storage.lookupRef<Cfg>(position).template get<TestDescriptor::Var4>() ==
            &ref.template get<TestDescriptor::Var4>());
    }
  }

  storage.allocateBuckets();
  storage.allocateScratchPads();
}

TEST_CASE("Storage scratchpad sharing" * doctest::test_suite("memory")) {
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2}),
      // the sharing does not depend on the configuration; any one built in does
      initializer::EnumLayer(std::vector<ConfigId>{defaultConfig()}));

  const auto setup = [&](initializer::Storage& storage) {
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
    initializer::Storage storage;
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
    initializer::Storage storage;
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

TEST_CASE("Type numbering" * doctest::test_suite("memory")) {
  // a numbering of its own: the ones of the process have given out ids already
  TypeNumbering numbering;
  CHECK(numbering.size() == 0);
  CHECK(numbering.idOf(typeid(int)) == 0);
  CHECK(numbering.idOf(typeid(float)) == 1);
  CHECK(numbering.idOf(typeid(int)) == 0);
  CHECK(numbering.size() == 2);

  // the ids of the process cover the types the program names, each with an id of its own
  const auto id1 = TypeId<TypeDomain::Storage, TestDescriptor::Var1>;
  const auto id2 = TypeId<TypeDomain::Storage, TestDescriptor::Var2>;
  CHECK(id1 != id2);
  CHECK(id1 < typeCount(TypeDomain::Storage));
  CHECK(id2 < typeCount(TypeDomain::Storage));
}

TEST_CASE_TEMPLATE("Storage views" * doctest::test_suite("memory"), Cfg, SEISSOL_CONFIG_TYPES) {
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1}),
      initializer::EnumLayer(std::vector<ConfigId>{configIdOf<Cfg>()}));

  initializer::Storage storage;
  // only on the interior layer
  storage.add<TestDescriptor::Var1>(initializer::LayerMask(Ghost) | initializer::LayerMask(Copy),
                                    sizeof(void*),
                                    initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var4>(
      initializer::LayerMask(), sizeof(void*), initializer::AllocationMode::HostOnly);
  CHECK(storage.has<TestDescriptor::Var1>());
  CHECK_FALSE(storage.has<TestDescriptor::Var2>());

  storage.setLayerCount(colorMap);
  storage.fixate();
  for (auto& layer : storage.leaves()) {
    layer.setNumberOfCells(3);
  }
  storage.allocateVariables();

  for (auto& layer : storage.leaves()) {
    auto ref = layer.cellRef<Cfg>(2);
    CHECK(&ref.template get<TestDescriptor::Var4>() == &layer.var<TestDescriptor::Var4>(Cfg())[2]);

    // what the layer does not hold is a null pointer: a filtered variable, and a variable that the
    // storage does not have
    const bool interior = layer.getIdentifier().halo == HaloType::Interior;
    CHECK((layer.view().address<TestDescriptor::Var1>(2) != nullptr) == interior);
    CHECK(layer.view().address<TestDescriptor::Var2>(2) == nullptr);

    // a derived view puts single variables elsewhere, and leaves the layer as it is
    PerConfigArray<Cfg> elsewhere{};
    initializer::LayerView view = layer.view();
    view.relocate<TestDescriptor::Var4>(&elsewhere);
    const initializer::CellRef<Cfg> moved(view, 2);
    CHECK(&moved.template get<TestDescriptor::Var4>() == &elsewhere);
    CHECK(&ref.template get<TestDescriptor::Var4>() == &layer.var<TestDescriptor::Var4>(Cfg())[2]);
    if (interior) {
      CHECK(&moved.template get<TestDescriptor::Var1>() ==
            &ref.template get<TestDescriptor::Var1>());
    }
  }
}

} // namespace seissol::unit_test
