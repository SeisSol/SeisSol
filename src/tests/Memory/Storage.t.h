// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <Common/ConfigHelper.h>
#include <Initializer/BasicTypedefs.h>
#include <Memory/Tree/Colormap.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
namespace seissol::unit_test {

using namespace seissol;

struct TestDescriptor {
  struct Var1 : public initializer::Variable<int> {};
  struct Var2 : public initializer::Variable<float> {};
  struct Var3 : public initializer::Variable<double[123]> {};
  struct Bucket : public initializer::Bucket<float> {};
  struct Scratchpad : public initializer::Scratchpad<float> {};
};

TEST_CASE("Storage") {
  initializer::Storage<initializer::GenericVarmap> storage;

  // NOTE: the LTSColorMap is hard-coded to the storage right now.
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1, 2, 3}),
      initializer::TraitLayer(std::vector<ConfigVariant>{ConfigVariantList[0]}));

  storage.add<TestDescriptor::Var1>(Ghost, 1, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var2>(Copy, 4, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Var3>(Interior, 8, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Bucket>(
      initializer::LayerMask(), 1, initializer::AllocationMode::HostOnly);
  storage.add<TestDescriptor::Scratchpad>(
      initializer::LayerMask(), 1, initializer::AllocationMode::HostOnly);

  storage.setLayerCount(colorMap);

  REQUIRE(storage.numChildren() == colorMap.size());

  storage.fixate();

  for (const auto [i, layer] : common::enumerate(storage.leaves())) {
    REQUIRE(layer.getIdentifier().lts == colorMap.argument(i).lts);
    REQUIRE(layer.getIdentifier().halo == colorMap.argument(i).halo);
    REQUIRE(layer.getIdentifier().config.index() == colorMap.argument(i).config.index());
  }

  for (auto [i, layer] : common::enumerate(storage.leaves())) {
    layer.setNumberOfCells(i + 1);
  }

  storage.allocateVariables();
  storage.touchVariables();

  for (const auto [i, layer] : common::enumerate(storage.leaves())) {
    REQUIRE(layer.size() == i + 1);
  }

  storage.allocateBuckets();
  storage.allocateScratchPads();
}

} // namespace seissol::unit_test
