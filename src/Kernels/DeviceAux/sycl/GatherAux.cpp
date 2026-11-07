// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/DeviceAux/GatherAux.h"

#include "Parallel/BlockSet.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <sycl/sycl.hpp>

namespace seissol::kernels::device::aux::gather {

namespace {
constexpr std::size_t Lanes = 32;
constexpr std::size_t MaxSets = 32;

using parallel::BlockSet;

struct Sets {
  BlockSet sets[MaxSets];
};
} // namespace

void copyBlocks(const parallel::BlockSet* sets,
                std::size_t setCount,
                void* const* blocks,
                void* packed,
                bool scatter,
                void* stream) {
  auto* queue = reinterpret_cast<sycl::queue*>(stream);
  for (std::size_t first = 0; first < setCount; first += MaxSets) {
    const auto count = std::min(MaxSets, setCount - first);
    Sets part{};
    std::size_t maxBlocks = 0;
    for (std::size_t set = 0; set < count; ++set) {
      part.sets[set] = sets[first + set];
      maxBlocks = std::max(maxBlocks, part.sets[set].count);
    }
    if (maxBlocks == 0) {
      continue;
    }
    auto* packedBytes = static_cast<unsigned char*>(packed);
    // the first dimension goes over the sets, the second over the lanes of their blocks; the lanes
    // of a block copy it in words of 4 bytes where both sides are aligned to them, in bytes
    // otherwise
    queue->parallel_for(sycl::range<2>(count, maxBlocks * Lanes), [=](sycl::item<2> item) {
      const BlockSet& set = part.sets[item[0]];
      const std::size_t i = item[1] / Lanes;
      const std::size_t lane = item[1] % Lanes;
      if (i >= set.count || blocks[set.begin + i] == nullptr) {
        return;
      }
      auto* remote = static_cast<unsigned char*>(blocks[set.begin + i]);
      unsigned char* local = packedBytes + set.offset + i * set.size;
      unsigned char* target = scatter ? remote : local;
      const unsigned char* source = scatter ? local : remote;
      const bool aligned = ((reinterpret_cast<std::uintptr_t>(target) |
                             reinterpret_cast<std::uintptr_t>(source) | set.size) %
                            sizeof(std::uint32_t)) == 0;
      if (aligned) {
        auto* targetWords = reinterpret_cast<std::uint32_t*>(target);
        const auto* sourceWords = reinterpret_cast<const std::uint32_t*>(source);
        for (std::size_t word = lane; word < set.size / sizeof(std::uint32_t); word += Lanes) {
          targetWords[word] = sourceWords[word];
        }
      } else {
        for (std::size_t byte = lane; byte < set.size; byte += Lanes) {
          target[byte] = source[byte];
        }
      }
    });
  }
}

} // namespace seissol::kernels::device::aux::gather
