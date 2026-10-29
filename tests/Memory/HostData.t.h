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

#include <array>
#include <cstddef>
#include <set>
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

  // the memory of each variable once
  const auto containers = some.containers();
  CHECK(containers.size() == 3);
  CHECK(std::set<const initializer::DualMemoryContainer*>(containers.begin(), containers.end())
            .size() == 3);
  // the copies of host-only variables do nothing
  some.copyToHost(nullptr);

  // all variables of a storage hold the single ones of it
  initializer::HostData all;
  all.add<TestDescriptor::Var1>(storage).addVariables(storage);
  CHECK(all.containers().size() == 2);

  some.merge(initializer::HostData::everything());
  CHECK(some.isEverything());
  CHECK_FALSE(some.empty());
  CHECK(initializer::HostData::everything().isEverything());
}

TEST_CASE("The snapshots of memory containers" * doctest::test_suite("memory")) {
  // only the pointers count here; nothing gets allocated or copied
  int hostPart = 0;
  int devicePart = 0;
  int snapshotPart = 0;

  initializer::DualMemoryContainer split;
  split.allocationMode = initializer::AllocationMode::HostDeviceSplitPinned;
  split.host = &hostPart;
  split.device = &devicePart;
  CHECK(split.separateParts());
  CHECK(split.copiesToHost());
  CHECK_FALSE(split.snapshotApartFromHost());
  // the snapshot on the device gets copied to the host part later
  split.snapshot = &snapshotPart;
  CHECK(split.snapshotApartFromHost());
  split.constant = true;
  CHECK_FALSE(split.copiesToHost());

  initializer::DualMemoryContainer unified;
  unified.allocationMode = initializer::AllocationMode::HostDeviceUnified;
  unified.host = &devicePart;
  unified.device = &devicePart;
  CHECK_FALSE(unified.separateParts());
  CHECK_FALSE(unified.copiesToHost());
  // as a host part of its own, the snapshot splits unified memory into two parts
  unified.host = &snapshotPart;
  unified.snapshot = &snapshotPart;
  CHECK(unified.separateParts());
  CHECK(unified.copiesToHost());
  CHECK_FALSE(unified.snapshotApartFromHost());

  // the layers see the snapshot at their offset, as they see the other parts
  std::array<char, 16> host{};
  std::array<char, 16> device{};
  std::array<char, 16> snapshot{};
  initializer::DualMemoryContainer whole;
  whole.allocationMode = initializer::AllocationMode::HostDeviceSplitPinned;
  whole.host = host.data();
  whole.device = device.data();
  whole.snapshot = snapshot.data();
  initializer::DualMemoryContainer part;
  part.offsetFrom(whole, 4, 8);
  CHECK(part.host == host.data() + 4);
  CHECK(part.device == device.data() + 4);
  CHECK(part.snapshot == snapshot.data() + 4);
  CHECK(part.allocationSize == 8);
}

} // namespace seissol::unit_test
