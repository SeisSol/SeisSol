// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Memory/HostSnapshot.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::unit_test {

#ifdef ACL_DEVICE

namespace hostsnapshot {

constexpr std::size_t Count = 1024;

/// writes `value` into all of the device part, as the device does while the time steps go on
inline void writeOnDevice(const initializer::DualMemoryContainer& container, int value) {
  const std::vector<int> values(Count, value);
  ::device::DeviceInstance::instance().api().copyTo(
      container.device, values.data(), Count * sizeof(int));
}

inline bool hostHolds(const initializer::DualMemoryContainer& container, int value) {
  const auto* host = static_cast<const int*>(container.host);
  return std::all_of(host, host + Count, [value](int entry) { return entry == value; });
}

} // namespace hostsnapshot

TEST_CASE("A snapshot gives the host the data of the time it was taken" *
          doctest::test_suite("memory")) {
  using namespace hostsnapshot;
  auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();
  void* stream = device.api().createStream();
  memory::ManagedAllocator allocator;
  constexpr std::size_t Alignment = 64;

  // pinned host memory: the snapshot is on the device, and gets copied to the host later
  initializer::DualMemoryContainer pinned;
  pinned.allocate(allocator,
                  Count * sizeof(int),
                  Alignment,
                  initializer::AllocationMode::HostDeviceSplitPinned);
  pinned.allocateSnapshot(allocator, Alignment);
  REQUIRE(pinned.snapshotApartFromHost());
  // host memory that is not pinned: no snapshot, the data gets to the host right away
  initializer::DualMemoryContainer split;
  split.allocate(
      allocator, Count * sizeof(int), Alignment, initializer::AllocationMode::HostDeviceSplit);
  split.allocateSnapshot(allocator, Alignment);
  REQUIRE(split.snapshot == nullptr);
  // unified memory: the snapshot is a host part of its own, apart from the device part
  initializer::DualMemoryContainer unified;
  unified.allocate(
      allocator, Count * sizeof(int), Alignment, initializer::AllocationMode::HostDeviceUnified);
  unified.allocateSnapshot(allocator, Alignment);
  REQUIRE(unified.snapshot == unified.host);
  REQUIRE(unified.host != unified.device);
  // constant data never gets copied back
  initializer::DualMemoryContainer constant;
  constant.constant = true;
  constant.allocate(allocator,
                    Count * sizeof(int),
                    Alignment,
                    initializer::AllocationMode::HostDeviceSplitPinned);
  std::fill_n(static_cast<int*>(constant.host), Count, 7);

  const std::vector<const initializer::DualMemoryContainer*> containers{
      &pinned, &split, &unified, &constant};
  for (const auto* container : containers) {
    writeOnDevice(*container, 1);
  }
  initializer::HostSnapshot snapshot;
  snapshot.take(containers, stream);
  // the device goes on, while the host has yet to read the data of the synchronization point
  for (const auto* container : containers) {
    writeOnDevice(*container, 2);
  }
  CHECK(hostHolds(split, 1));
  CHECK(hostHolds(unified, 1));

  // e.g. on the thread of the outputs
  std::thread([&snapshot]() { snapshot.copyToHost(); }).join();
  CHECK(hostHolds(pinned, 1));
  CHECK(hostHolds(split, 1));
  CHECK(hostHolds(unified, 1));
  CHECK(hostHolds(constant, 7));

  // the next synchronization point
  snapshot.take(containers, stream);
  snapshot.copyToHost();
  CHECK(hostHolds(pinned, 2));
  CHECK(hostHolds(split, 2));
  CHECK(hostHolds(unified, 2));
  CHECK(hostHolds(constant, 7));

  device.api().destroyGenericStream(stream);
}

#endif

} // namespace seissol::unit_test
