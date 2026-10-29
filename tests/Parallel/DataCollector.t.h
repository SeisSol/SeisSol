// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Kernels/Common.h"
#include "Memory/MemoryAllocator.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Runtime/HostBufferPool.h"
#include "Parallel/Runtime/Stream.h"

#include <cstddef>
#include <cstdint>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::unit_test {

namespace datacollectortest {
inline void setUpDevice() {
#ifdef ACL_DEVICE
  auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();
#endif
}

// the memory that the blocks live in: on the device, or (without one) where the host reads them
constexpr auto BlockMemkind =
    isDeviceOn() ? memory::Memkind::DeviceGlobalMemory : memory::Memkind::Standard;

inline std::uint8_t
    pattern(std::size_t set, std::size_t block, std::size_t byte, std::size_t round) {
  return static_cast<std::uint8_t>((set * 31 + block * 7 + byte + round * 101) % 251);
}

/**
 * Sets of blocks in one allocation per set, at consecutive offsets (so that blocks of sizes that
 * are no multiple of a word are not aligned to one either); every third block of a set with more
 * than two blocks is left out (null).
 */
struct Blocks {
  std::vector<std::size_t> sizes;
  std::vector<std::vector<void*>> sets;
  std::vector<void*> allocations;

  explicit Blocks(const std::vector<std::size_t>& blockSizes, std::size_t count)
      : sizes(blockSizes) {
    for (std::size_t set = 0; set < sizes.size(); ++set) {
      auto* memory =
          static_cast<std::uint8_t*>(memory::allocate(count * sizes[set], 64, BlockMemkind));
      allocations.push_back(memory);
      std::vector<std::uint8_t> host(count * sizes[set]);
      std::vector<void*> blocks;
      for (std::size_t block = 0; block < count; ++block) {
        for (std::size_t byte = 0; byte < sizes[set]; ++byte) {
          host[block * sizes[set] + byte] = pattern(set, block, byte, 0);
        }
        const bool skipped = count > 2 && block % 3 == 1;
        blocks.push_back(skipped ? nullptr : memory + block * sizes[set]);
      }
      memory::memcopy(memory, host.data(), host.size(), BlockMemkind, memory::Memkind::Standard);
      sets.push_back(blocks);
    }
  }

  ~Blocks() {
    for (auto* allocation : allocations) {
      memory::free(allocation, BlockMemkind);
    }
  }

  Blocks(const Blocks&) = delete;
  Blocks(Blocks&&) = delete;
  Blocks& operator=(const Blocks&) = delete;
  Blocks& operator=(Blocks&&) = delete;

  /// whether the block holds the pattern of `round` (read back from where it lives)
  [[nodiscard]] bool holds(std::size_t set, std::size_t block, std::size_t round) const {
    std::vector<std::uint8_t> host(sizes[set]);
    memory::memcopy(
        host.data(), sets[set][block], sizes[set], memory::Memkind::Standard, BlockMemkind);
    for (std::size_t byte = 0; byte < sizes[set]; ++byte) {
      if (host[byte] != pattern(set, block, byte, round)) {
        return false;
      }
    }
    return true;
  }
};

inline bool holds(
    const void* data, std::size_t size, std::size_t set, std::size_t block, std::size_t round) {
  const auto* bytes = static_cast<const std::uint8_t*>(data);
  for (std::size_t byte = 0; byte < size; ++byte) {
    if (bytes[byte] != pattern(set, block, byte, round)) {
      return false;
    }
  }
  return true;
}
} // namespace datacollectortest

using namespace datacollectortest;

TEST_CASE("A data collector gathers and scatters sets of blocks of different sizes" *
          doctest::test_suite("parallel")) {
  setUpDevice();
  // more sets than one kernel takes, of sizes in words and of sizes in bytes only
  constexpr std::size_t SetCount = 40;
  std::vector<std::size_t> sizes;
  sizes.reserve(SetCount);
  for (std::size_t set = 0; set < SetCount; ++set) {
    sizes.push_back(set % 5 == 4 ? 6 + set : 4 * (1 + set % 7) * 16);
  }
  const Blocks blocks(sizes, 5);
  parallel::DataCollectorUntyped collector(blocks.sets, blocks.sizes);
  parallel::runtime::StreamRuntime runtime;

  collector.gatherToHost(runtime.stream());
  runtime.wait();
  for (std::size_t set = 0; set < sizes.size(); ++set) {
    for (std::size_t block = 0; block < blocks.sets[set].size(); ++block) {
      if (blocks.sets[set][block] != nullptr) {
        CHECK(holds(collector.get(set, block), sizes[set], set, block, 0));
      }
    }
  }

  // the host changes the gathered data, which goes back to the blocks
  for (std::size_t set = 0; set < sizes.size(); ++set) {
    for (std::size_t block = 0; block < blocks.sets[set].size(); ++block) {
      if (blocks.sets[set][block] != nullptr) {
        auto* bytes = static_cast<std::uint8_t*>(collector.get(set, block));
        for (std::size_t byte = 0; byte < sizes[set]; ++byte) {
          bytes[byte] = pattern(set, block, byte, 1);
        }
      }
    }
  }
  collector.scatterFromHost(runtime.stream());
  runtime.wait();
  for (std::size_t set = 0; set < sizes.size(); ++set) {
    for (std::size_t block = 0; block < blocks.sets[set].size(); ++block) {
      if (blocks.sets[set][block] != nullptr) {
        CHECK(blocks.holds(set, block, 1));
      }
    }
  }
  runtime.dispose();
}

TEST_CASE("A data collector gathers into buffers of a host buffer pool, one per round" *
          doctest::test_suite("parallel")) {
  setUpDevice();
  const Blocks blocks({2016, 6, 64}, 4);
  parallel::DataCollectorUntyped collector(blocks.sets, blocks.sizes);
  parallel::runtime::HostBufferPool pool(collector.gatheredBytes());
  parallel::runtime::StreamRuntime runtime;

  // two rounds in flight at once: each one keeps what the device gathered for it
  const auto first = pool.take();
  collector.gatherTo(first.device, runtime.stream());
  const auto second = pool.take();
  CHECK(second.host != first.host);
  collector.gatherTo(second.device, runtime.stream());
  runtime.wait();
  for (const auto& buffer : {first, second}) {
    for (std::size_t set = 0; set < blocks.sizes.size(); ++set) {
      for (std::size_t block = 0; block < blocks.sets[set].size(); ++block) {
        if (blocks.sets[set][block] != nullptr) {
          CHECK(holds(collector.get(set, block, buffer.host), blocks.sizes[set], set, block, 0));
        }
      }
    }
  }
  pool.give(first);
  pool.give(second);
  CHECK(pool.size() == 2);
  runtime.dispose();
}

TEST_CASE("A data collector of a single set gathers its blocks" * doctest::test_suite("parallel")) {
  setUpDevice();
  const Blocks blocks({2016}, 2);
  parallel::DataCollectorUntyped collector(blocks.sets[0], blocks.sizes[0]);
  parallel::runtime::StreamRuntime runtime;
  collector.gatherToHost(runtime.stream());
  runtime.wait();
  CHECK(holds(collector.get(0), blocks.sizes[0], 0, 0, 0));
  CHECK(holds(collector.get(1), blocks.sizes[0], 0, 1, 0));
  runtime.dispose();
}

} // namespace seissol::unit_test
