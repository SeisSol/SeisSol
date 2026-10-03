// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Iterator.h"
#include "Common/Real.h"
#include "Config.h"
#include "Initializer/BasicTypedefs.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"

#include <cstddef>
#include <type_traits>
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

TEST_CASE("Storage" * doctest::test_suite("memory")) {
  initializer::Storage<initializer::GenericVarmap> storage;

  // NOTE: the LTSColorMap is hard-coded to the storage right now.
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2, 3}),
      initializer::EnumLayer(std::vector<ConfigId>{configIdOf<Config>()}));

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
          sizeof(PerConfigArray<Config>));
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
    auto* perConfig = layer.var<TestDescriptor::Var4>(Config());
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t j = 0; j <= Config::ConvergenceOrder; ++j) {
        CHECK(perConfig[cell][j] == 0);
        perConfig[cell][j] = static_cast<RealT<Config::Precision>>(cell + j);
      }
    }
  }

  // a cell of a configuration sees the variables as that configuration holds them, whichever way
  // it is reached
  for (auto [color, layer] : common::enumerate(storage.leaves())) {
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      auto ref = layer.cellRef<Config>(cell);
      static_assert(
          std::is_same_v<decltype(ref.get<TestDescriptor::Var4>()), PerConfigArray<Config>&>);
      static_assert(std::is_same_v<decltype(ref.get<TestDescriptor::Var3>()), double (&)[123]>);
      CHECK(&ref.get<TestDescriptor::Var4>() == &layer.var<TestDescriptor::Var4>(Config())[cell]);

      initializer::StoragePosition position;
      position.color = color;
      position.cell = cell;
      CHECK(&storage.lookup<TestDescriptor::Var4>(Config(), position) ==
            &ref.get<TestDescriptor::Var4>());
      CHECK(&storage.lookupRef<Config>(position).get<TestDescriptor::Var4>() ==
            &ref.get<TestDescriptor::Var4>());
    }
  }

  storage.allocateBuckets();
  storage.allocateScratchPads();
}

} // namespace seissol::unit_test
