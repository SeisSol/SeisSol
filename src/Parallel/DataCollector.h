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

#include <cstddef>
#include <vector>

namespace seissol::parallel {

// wrapper class for sparse device-to-host transfers (e.g. needed for receivers)
class DataCollectorUntyped {
  public:
  DataCollectorUntyped(const std::vector<void*>& indexDataHost,
                       size_t elemSize,
                       bool hostAccessible = false);

  ~DataCollectorUntyped();

  auto operator=(const DataCollectorUntyped&) = delete;
  DataCollectorUntyped(const DataCollectorUntyped&) = delete;

  auto operator=(DataCollectorUntyped&&) -> DataCollectorUntyped& = default;
  DataCollectorUntyped(DataCollectorUntyped&&) = default;

  // gathers and sends data to the device to the host
  void gatherToHost(void* stream);

  // sends and scatters data from the host to the device
  void scatterFromHost(void* stream);

  void* get(size_t index) {
    if (hostAccessible_) {
      return indexDataHost_[index];
    } else {
      return reinterpret_cast<void*>(reinterpret_cast<uint8_t*>(copiedData_) + index * elemSize_);
    }
  }

  [[nodiscard]] const void* get(size_t index) const {
    if (hostAccessible_) {
      return indexDataHost_[index];
    } else {
      return reinterpret_cast<void*>(reinterpret_cast<uint8_t*>(copiedData_) + index * elemSize_);
    }
  }

  private:
  bool hostAccessible_;
  void** indexDataDevice_;
  std::vector<void*> indexDataHost_;
  size_t indexCount_;
  size_t elemSize_;
  void* copiedData_;
  void* copiedDataDevice_;
};

// wrapper class; typed version
template <typename T>
class DataCollector : DataCollectorUntyped {
  public:
  DataCollector(const std::vector<T*>& data, std::size_t elemCount, bool hostAccessible = false)
      : DataCollectorUntyped(
            std::vector<void*>(data.begin(), data.end()), elemCount * sizeof(T), hostAccessible) {}

  void gatherToHost(void* stream) { DataCollectorUntyped::gatherToHost(stream); }

  void scatterFromHost(void* stream) { DataCollectorUntyped::scatterFromHost(stream); }

  T* get(std::size_t index) { return reinterpret_cast<T*>(DataCollectorUntyped::get(index)); }

  [[nodiscard]] const T* get(std::size_t index) const {
    return reinterpret_cast<const T*>(DataCollectorUntyped::get(index));
  }
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_DATACOLLECTOR_H_
