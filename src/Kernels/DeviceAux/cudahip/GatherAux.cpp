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

#ifdef __HIP__
#include "hip/hip_runtime.h"
using StreamT = hipStream_t;
// HIP keeps the arguments of a kernel in memory that it can index anyway
#define SEISSOL_GRID_CONSTANT
#endif
#ifdef __CUDACC__
using StreamT = cudaStream_t;
// indexed by the block, the sets stay where the arguments are instead of getting copied per thread
#define SEISSOL_GRID_CONSTANT __grid_constant__
#endif

namespace seissol::kernels::device::aux::gather {

namespace {
constexpr unsigned Lanes = 32;
// the blocks that a group of threads copies at a time, one per warp
constexpr unsigned BlocksPerGroup = 8;
constexpr std::size_t MaxGroups = 256;
constexpr std::size_t MaxSets = 32;

using parallel::BlockSet;

struct Sets {
  BlockSet sets[MaxSets];
};

// the second grid dimension goes over the sets, the first over their blocks; a warp copies a block,
// in words of 4 bytes where both sides are aligned to them, and in bytes otherwise
template <bool Scatter>
__global__ void copyBlocksKernel(const SEISSOL_GRID_CONSTANT Sets sets,
                                 void* const* __restrict blocks,
                                 unsigned char* __restrict packed) {
  const BlockSet& set = sets.sets[blockIdx.y];
  const std::size_t stride = static_cast<std::size_t>(gridDim.x) * blockDim.y;
  for (std::size_t i = static_cast<std::size_t>(blockIdx.x) * blockDim.y + threadIdx.y;
       i < set.count;
       i += stride) {
    auto* remote = static_cast<unsigned char*>(blocks[set.begin + i]);
    if (remote == nullptr) {
      continue;
    }
    unsigned char* local = packed + set.offset + i * set.size;
    unsigned char* target = Scatter ? remote : local;
    const unsigned char* source = Scatter ? local : remote;
    const bool aligned = ((reinterpret_cast<std::uintptr_t>(target) |
                           reinterpret_cast<std::uintptr_t>(source) | set.size) %
                          sizeof(std::uint32_t)) == 0;
    if (aligned) {
      auto* targetWords = reinterpret_cast<std::uint32_t*>(target);
      const auto* sourceWords = reinterpret_cast<const std::uint32_t*>(source);
      for (std::size_t word = threadIdx.x; word < set.size / sizeof(std::uint32_t);
           word += blockDim.x) {
        targetWords[word] = sourceWords[word];
      }
    } else {
      for (std::size_t byte = threadIdx.x; byte < set.size; byte += blockDim.x) {
        target[byte] = source[byte];
      }
    }
  }
}
} // namespace

void copyBlocks(const parallel::BlockSet* sets,
                std::size_t setCount,
                void* const* blocks,
                void* packed,
                bool scatter,
                void* stream) {
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
    const dim3 threads(Lanes, BlocksPerGroup);
    const dim3 grid(
        std::clamp<std::size_t>((maxBlocks + BlocksPerGroup - 1) / BlocksPerGroup, 1, MaxGroups),
        count);
    auto* packedBytes = static_cast<unsigned char*>(packed);
    if (scatter) {
      copyBlocksKernel<true>
          <<<grid, threads, 0, reinterpret_cast<StreamT>(stream)>>>(part, blocks, packedBytes);
    } else {
      copyBlocksKernel<false>
          <<<grid, threads, 0, reinterpret_cast<StreamT>(stream)>>>(part, blocks, packedBytes);
    }
  }
}

} // namespace seissol::kernels::device::aux::gather
