// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_
#define SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_

#include "Kernels/Common.h"
#include "Memory/MemoryAllocator.h"
#include "Parallel/BlockSet.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::parallel {

/**
 * Sparse transfers from the device to the host (e.g. needed for receivers): the data at the given
 * device addresses gets gathered, in the order of the stream of whoever changes it, into pinned
 * memory on the host, by a kernel. If the host can read the memory of the device
 * (`hostAccessible`, e.g. unified memory), it gets gathered into a buffer of that memory instead,
 * which the host reads, so that it never reads what the device works on. Without a device, the
 * host reads the data where it is.
 */
class DataCollectorUntyped {
  public:
  DataCollectorUntyped(const std::vector<void*>& indexDataHost,
                       size_t elemSize,
                       bool hostAccessible = false);

  /**
   * Gathers several sets of blocks with one kernel (e.g. all variables of the fault output): the
   * blocks of set `s` are at `sets[s]`, with `sizes[s]` bytes each; get(s, index) gives block
   * `index` of set `s`.
   */
  DataCollectorUntyped(const std::vector<std::vector<void*>>& sets,
                       const std::vector<size_t>& sizes,
                       bool hostAccessible = false);

  ~DataCollectorUntyped();

  auto operator=(const DataCollectorUntyped&) = delete;
  DataCollectorUntyped(const DataCollectorUntyped&) = delete;

  auto operator=(DataCollectorUntyped&&) -> DataCollectorUntyped& = default;
  DataCollectorUntyped(DataCollectorUntyped&&) = default;

  // gathers and sends data to the device to the host
  void gatherToHost(void* stream);

  /**
   * Gathers the data into `destination` instead (as the device addresses it), which has room for
   * gatheredBytes() (e.g. a buffer of a runtime::HostBufferPool); get(set, index, gathered) finds
   * the blocks there. Without a device, there is nothing to gather.
   */
  void gatherTo(void* destination, void* stream);

  /// the bytes of all blocks, as they get gathered
  [[nodiscard]] size_t gatheredBytes() const { return gatheredBytes_; }

  /**
   * Block `index` of set `set` in what gatherTo() has gathered into `gathered` (as the host
   * addresses it); without a device, the block where it is.
   */
  [[nodiscard]] const void* get(size_t set, size_t index, const void* gathered) const;

  // sends and scatters data from the host to the device
  void scatterFromHost(void* stream);

  void* get(size_t index) { return get(0, index); }

  [[nodiscard]] const void* get(size_t index) const { return get(0, index); }

  void* get(size_t set, size_t index) {
    const auto& blocks = sets_[set];
    if (direct_) {
      return indexDataHost_[blocks.begin + index];
    }
    return static_cast<uint8_t*>(copiedData_) + blocks.offset + index * blocks.size;
  }

  [[nodiscard]] const void* get(size_t set, size_t index) const {
    const auto& blocks = sets_[set];
    if (direct_) {
      return indexDataHost_[blocks.begin + index];
    }
    return static_cast<const uint8_t*>(copiedData_) + blocks.offset + index * blocks.size;
  }

  private:
  // without a device: the host reads the data where it is
  bool direct_{false};
  bool hostAccessible_;
  void** indexDataDevice_{nullptr};
  // the blocks of all sets, one after the other
  std::vector<void*> indexDataHost_;
  size_t indexCount_{0};
  // the sets, and where their blocks lie in the gathered data
  std::vector<BlockSet> sets_;
  size_t gatheredBytes_{0};
  // the gathered data: pinned memory on the host, or a buffer of memory that both can read
  void* copiedData_{nullptr};
  void* copiedDataDevice_{nullptr};
};

// wrapper class; typed version
template <typename T>
class DataCollector : DataCollectorUntyped {
  public:
  DataCollector(const std::vector<T*>& data, std::size_t elemCount, bool hostAccessible = false)
      : DataCollectorUntyped(
            std::vector<void*>(data.begin(), data.end()), elemCount * sizeof(T), hostAccessible) {}

  void gatherToHost(void* stream) { DataCollectorUntyped::gatherToHost(stream); }

  void gatherTo(void* destination, void* stream) {
    DataCollectorUntyped::gatherTo(destination, stream);
  }

  [[nodiscard]] size_t gatheredBytes() const { return DataCollectorUntyped::gatheredBytes(); }

  [[nodiscard]] const T* get(std::size_t index, const void* gathered) const {
    return reinterpret_cast<const T*>(DataCollectorUntyped::get(0, index, gathered));
  }

  void scatterFromHost(void* stream) { DataCollectorUntyped::scatterFromHost(stream); }

  T* get(std::size_t index) { return reinterpret_cast<T*>(DataCollectorUntyped::get(index)); }

  [[nodiscard]] const T* get(std::size_t index) const {
    return reinterpret_cast<const T*>(DataCollectorUntyped::get(index));
  }
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_
