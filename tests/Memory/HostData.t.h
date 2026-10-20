// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigRegistry.h"
#include "Initializer/BasicTypedefs.h"
#include "Memory/HostData.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"

#include <cstddef>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("HostData" * doctest::test_suite("memory")) {
  const initializer::LTSColorMap colorMap(
      initializer::EnumLayer(
          std::vector<HaloType>{HaloType::Interior, HaloType::Copy, HaloType::Ghost}),
      initializer::EnumLayer(std::vector<std::size_t>{1}),
      initializer::EnumLayer(std::vector<ConfigId>{defaultConfig()}));

  const auto setup = [&](initializer::Storage<initializer::GenericVarmap>& storage) {
    storage.add<TestDescriptor::Var1>(
        initializer::LayerMask(), sizeof(void*), initializer::AllocationMode::HostOnly);
    storage.add<TestDescriptor::Var2>(
        initializer::LayerMask(), sizeof(void*), initializer::AllocationMode::HostOnly);
    storage.setLayerCount(colorMap);
    storage.fixate();
    for (auto& layer : storage.leaves()) {
      layer.setNumberOfCells(2);
    }
    storage.allocateVariables();
  };
  initializer::Storage<initializer::GenericVarmap> storage;
  initializer::Storage<initializer::GenericVarmap> otherStorage;
  setup(storage);
  setup(otherStorage);

  const initializer::HostData none;
  CHECK(none.empty());
  CHECK_FALSE(none.isEverything());
  CHECK(none.size() == 0);

  // a variable counts once, however often it is asked for; the same variable of another storage
  // counts separately
  initializer::HostData some;
  some.add<TestDescriptor::Var1>(storage)
      .add<TestDescriptor::Var1>(storage)
      .add<TestDescriptor::Var2>(storage);
  CHECK(some.size() == 2);
  CHECK_FALSE(some.empty());

  initializer::HostData other;
  other.add<TestDescriptor::Var2>(storage).add<TestDescriptor::Var2>(otherStorage);
  some.merge(other);
  CHECK(some.size() == 3);
  CHECK_FALSE(some.isEverything());

  // host-only variables need no copies
  CHECK(some.copies().empty());
  some.copyToHost(nullptr);

  some.merge(initializer::HostData::everything());
  CHECK(some.isEverything());
  CHECK_FALSE(some.empty());
  CHECK(initializer::HostData::everything().isEverything());
}

} // namespace seissol::unit_test
