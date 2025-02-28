// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LAYER_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LAYER_H_

#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/DeviceGraph.h"
#include "Memory/MemoryAllocator.h"
#include "Node.h"
#include <bitset>
#include <cstring>
#include <limits>
#include <type_traits>
#include <vector>

enum LayerType { Ghost = (1 << 0), Copy = (1 << 1), Interior = (1 << 2), NumLayers = 3 };

namespace seissol::initializer {
using LayerMask = std::bitset<NumLayers>;

enum class AllocationMode {
  HostOnly,
  HostOnlyHBM,
  HostDeviceUnified,
  HostDeviceSplit,
  HostDeviceSplitPinned,
  HostDevicePinned,
  DeviceOnly
};

enum class AllocationPlace { Host, Device };

struct DualMemoryContainer {
  void* host = nullptr;
  void* device = nullptr;
  AllocationMode allocationMode;
  std::size_t allocationSize{};
  std::size_t allocationAlignment{};
  bool constant{};

  [[nodiscard]] void* get(AllocationPlace place) const {
    if (place == AllocationPlace::Host) {
      return host;
    } else if (place == AllocationPlace::Device) {
      return device;
    } else {
      return nullptr; // should not happen
    }
  }

  void allocate(seissol::memory::ManagedAllocator& allocator,
                std::size_t size,
                std::size_t alignment,
                AllocationMode mode) {
    if (mode == AllocationMode::HostOnly) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::Standard);
      device = nullptr;
    }
    if (mode == AllocationMode::HostOnlyHBM) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::HighBandwidth);
      device = nullptr;
    }
    if (mode == AllocationMode::DeviceOnly) {
      host = nullptr;
      device =
          allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceGlobalMemory);
    }
    if (mode == AllocationMode::HostDeviceUnified) {
      host =
          allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceUnifiedMemory);
      device = host;
    }
    if (mode == AllocationMode::HostDevicePinned) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::PinnedMemory);
      device = seissol::memory::hostToDevicePointer(host, seissol::memory::Memkind::PinnedMemory);
    }
    if (mode == AllocationMode::HostDeviceSplit) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::Standard);
      device =
          allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceGlobalMemory);
    }
    if (mode == AllocationMode::HostDeviceSplitPinned) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::PinnedMemory);
      device =
          allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceGlobalMemory);
    }
    allocationMode = mode;
    allocationSize = size;
  }

  void synchronizeTo(AllocationPlace place, void* stream) {
#ifdef ACL_DEVICE
    if (allocationMode == AllocationMode::HostDeviceSplit ||
        allocationMode == AllocationMode::HostDeviceSplitPinned) {
      if (place == AllocationPlace::Host) {
        // do not copy back constant data (we ignore the other direction for now)
        if (!constant) {
          device::DeviceInstance::getInstance().api->copyFromAsync(
              host, device, allocationSize, stream);
        }
      } else {
        device::DeviceInstance::getInstance().api->copyToAsync(
            device, host, allocationSize, stream);
      }
    }
    if (allocationMode == AllocationMode::HostDeviceUnified) {
      // currently broken (?)
      if (place == AllocationPlace::Host) {
        // device::DeviceInstance::getInstance().api->prefetchUnifiedMemTo(device::Destination::Host,
        // host, allocationSize, stream);
      } else {
        // device::DeviceInstance::getInstance().api->prefetchUnifiedMemTo(device::Destination::CurrentDevice,
        // device, allocationSize, stream);
      }
    }
#endif
  }

  void offsetFrom(const DualMemoryContainer& container, std::size_t offset, std::size_t size) {
    host = nullptr;
    device = nullptr;

    allocationMode = container.allocationMode;
    allocationSize = size;
    if (container.host != nullptr) {
      host = reinterpret_cast<char*>(container.host) + offset;
    }
    if (container.device != nullptr) {
      device = reinterpret_cast<char*>(container.device) + offset;
    }
  }
};

template <typename T>
struct Variable {
  unsigned index;
  LayerMask mask;
  unsigned count{1};
  Variable() : index(std::numeric_limits<unsigned>::max()) {}
};

struct Bucket {
  unsigned index;

  Bucket() : index(std::numeric_limits<unsigned>::max()) {}
};

#ifdef ACL_DEVICE
struct ScratchpadMemory : public Bucket {};
#endif

struct MemoryInfo {
  size_t bytes{};
  size_t alignment{};
  size_t elemsize{};
  LayerMask mask;
  // seissol::memory::Memkind memkind;
  AllocationMode allocMode;
  bool constant{false};
};

class Layer : public Node {
  private:
  enum LayerType m_layerType;
  unsigned m_numberOfCells{0};
  std::vector<DualMemoryContainer> m_vars;
  std::vector<DualMemoryContainer> m_buckets;
  std::vector<size_t> m_bucketSizes;

#ifdef ACL_DEVICE
  std::vector<DualMemoryContainer> m_scratchpads{};
  std::vector<size_t> m_scratchpadSizes{};
  std::unordered_map<GraphKey, device::DeviceGraphHandle, GraphKeyHash> m_computeGraphHandles{};
  ConditionalPointersToRealsTable m_conditionalPointersToRealsTable{};
  DrConditionalPointersToRealsTable m_drConditionalPointersToRealsTable{};
  ConditionalMaterialTable m_conditionalMaterialTable{};
  ConditionalIndicesTable m_conditionalIndicesTable;
#endif

  public:
  Layer() = default;
  ~Layer() override = default;

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& variable : m_vars) {
      variable.synchronizeTo(place, stream);
    }
    for (auto& bucket : m_buckets) {
      bucket.synchronizeTo(place, stream);
    }
#ifdef ACL_DEVICE
    for (auto& scratchpad : m_scratchpads) {
      scratchpad.synchronizeTo(place, stream);
    }
#endif
  }

  template <typename T>
  T* var(const Variable<T>& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_vars.size() > handle.index);
    return static_cast<T*>(m_vars[handle.index].get(place));
  }

  template <typename T>
  void varSynchronizeTo(const Variable<T>& handle, AllocationPlace place, void* stream) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_vars.size() > handle.index);
    m_vars[handle.index].synchronizeTo(place, stream);
  }

  void* bucket(const Bucket& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_buckets.size() > handle.index);
    return m_buckets[handle.index].get(place);
  }

  void bucketSynchronizeTo(const Bucket& handle, AllocationPlace place, void* stream) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_buckets.size() > handle.index);
    m_buckets[handle.index].synchronizeTo(place, stream);
  }

#ifdef ACL_DEVICE
  void* getScratchpadMemory(const ScratchpadMemory& handle,
                            AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_scratchpads.size() > handle.index);
    return (m_scratchpads[handle.index].get(place));
  }
#endif

  /// i-th bit of layerMask shall be set if data is masked on the i-th layer
  [[nodiscard]] bool isMasked(LayerMask layerMask) const {
    return (LayerMask(m_layerType) & layerMask).any();
  }

  void setLayerType(enum LayerType layerType) { m_layerType = layerType; }

  [[nodiscard]] enum LayerType getLayerType() const { return m_layerType; }

  [[nodiscard]] unsigned getNumberOfCells() const { return m_numberOfCells; }

  void setNumberOfCells(unsigned numberOfCells) { m_numberOfCells = numberOfCells; }

  void allocatePointerArrays(unsigned numVars, unsigned numBuckets) {
    assert(m_vars.empty() && m_buckets.empty() && m_bucketSizes.empty());
    m_vars.resize(numVars);
    m_buckets.resize(numBuckets);
    m_bucketSizes.resize(numBuckets, 0);
  }

#ifdef ACL_DEVICE
  inline void allocateScratchpadArrays(unsigned numScratchPads) {
    assert(m_scratchpads.empty() && m_scratchpadSizes.empty());

    m_scratchpads.resize(numScratchPads);
    m_scratchpadSizes.resize(numScratchPads, 0);
  }
#endif

  void setBucketSize(const Bucket& handle, size_t size) {
    assert(m_bucketSizes.size() > handle.index);
    m_bucketSizes[handle.index] = size;
  }

#ifdef ACL_DEVICE
  inline void setScratchpadSize(const ScratchpadMemory& handle, size_t size) {
    assert(m_scratchpadSizes.size() > handle.index);
    m_scratchpadSizes[handle.index] = size;
  }
#endif

  size_t getBucketSize(const Bucket& handle) {
    assert(m_bucketSizes.size() > handle.index);
    return m_bucketSizes[handle.index];
  }

  void addVariableSizes(const std::vector<MemoryInfo>& vars, std::vector<size_t>& bytes) const {
    for (unsigned var = 0; var < vars.size(); ++var) {
      if (!isMasked(vars[var].mask)) {
        bytes[var] += m_numberOfCells * vars[var].bytes;
      }
    }
  }

  void addBucketSizes(std::vector<size_t>& bytes) {
    for (unsigned bucket = 0; bucket < bytes.size(); ++bucket) {
      bytes[bucket] += m_bucketSizes[bucket];
    }
  }

#ifdef ACL_DEVICE
  // Overrides array's elements; if the corresponding local
  // scratchpad mem. size is bigger then the one inside of the array
  void findMaxScratchpadSizes(std::vector<size_t>& bytes) {
    for (size_t id = 0; id < bytes.size(); ++id) {
      bytes[id] = std::max(bytes[id], m_scratchpadSizes[id]);
    }
  }
#endif

  void setMemoryRegionsForVariables(const std::vector<MemoryInfo>& vars,
                                    const std::vector<DualMemoryContainer>& memory,
                                    const std::vector<size_t>& offsets) {
    assert(m_vars.size() >= vars.size());
    for (unsigned var = 0; var < vars.size(); ++var) {
      if (!isMasked(vars[var].mask)) {
        m_vars[var].offsetFrom(memory[var], offsets[var], m_numberOfCells * vars[var].bytes);
      }
    }
  }

  void setMemoryRegionsForBuckets(const std::vector<DualMemoryContainer>& memory,
                                  const std::vector<size_t>& offsets) {
    assert(m_buckets.size() >= offsets.size());
    for (unsigned bucket = 0; bucket < offsets.size(); ++bucket) {
      m_buckets[bucket].offsetFrom(memory[bucket], offsets[bucket], m_bucketSizes[bucket]);
    }
  }

#ifdef ACL_DEVICE
  void setMemoryRegionsForScratchpads(const std::vector<DualMemoryContainer>& memory) {
    assert(m_scratchpads.size() == memory.size());
    for (size_t id = 0; id < m_scratchpads.size(); ++id) {
      m_scratchpads[id] = memory[id];
    }
  }
#endif

  void touchVariables(const std::vector<MemoryInfo>& vars) {
    for (unsigned var = 0; var < vars.size(); ++var) {

      // NOTE: we don't touch device global memory because it is in a different address space
      // we will do deep-copy from the host to a device later on
      if (!isMasked(vars[var].mask) && (m_vars[var].host != nullptr)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (unsigned cell = 0; cell < m_numberOfCells; ++cell) {
          memset(static_cast<char*>(m_vars[var].host) + cell * vars[var].bytes, 0, vars[var].bytes);
        }
      }
    }
  }

#ifdef ACL_DEVICE
  template <typename InnerKeyType>
  auto& getConditionalTable() {
    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Wp>) {
      return m_conditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Dr>) {
      return m_drConditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Material>) {
      return m_conditionalMaterialTable;
    }

    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Indices>) {
      return m_conditionalIndicesTable;
    }
  }

  template <typename InnerKeyType>
  const auto& getConditionalTable() const {
    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Wp>) {
      return m_conditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Dr>) {
      return m_drConditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Material>) {
      return m_conditionalMaterialTable;
    }

    if constexpr (std::is_same_v<InnerKeyType, inner_keys::Indices>) {
      return m_conditionalIndicesTable;
    }
  }

  device::DeviceGraphHandle getDeviceComputeGraphHandle(GraphKey graphKey) {
    if (m_computeGraphHandles.find(graphKey) != m_computeGraphHandles.end()) {
      return m_computeGraphHandles[graphKey];
    } else {
      return device::DeviceGraphHandle();
    }
  }

  void updateDeviceComputeGraphHandle(GraphKey graphKey, device::DeviceGraphHandle graphHandle) {
    assert(m_computeGraphHandles.find(graphKey) == m_computeGraphHandles.end() &&
           "an entry of hash table must be empty on write");
    if (graphHandle.isInitialized()) {
      m_computeGraphHandles[graphKey] = graphHandle;
    }
  }
#endif // ACL_DEVICE
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TREE_LAYER_H_
