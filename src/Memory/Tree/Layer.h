// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_TREE_LAYER_H_
#define SEISSOL_SRC_MEMORY_TREE_LAYER_H_

#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/DeviceGraph.h"
#include "Memory/MemoryAllocator.h"
#include "Node.h"
#include <bitset>
#include <cstring>
#include <limits>
#include <type_traits>
#include <typeindex>
#include <variant>
#include <vector>

#include <utils/logger.h>

enum LayerType { Ghost = (1U << 0U), Copy = (1U << 1U), Interior = (1U << 2U), NumLayers = 3U };

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

struct ScratchpadMemory : public Bucket {};

enum class MemoryType { Variable, Bucket, Scratchpad };

struct MemoryInfo {
  size_t bytes{};
  size_t alignment{};
  size_t size{};
  LayerMask mask;
  AllocationMode allocMode;
  MemoryType type;
  bool constant{false};
  bool filtered{false};
};

struct LayerIdentifier {
  enum LayerType halo;
  int lts;
  std::variant<Config> config;
};

class Layer : public Node {
  private:
  LayerIdentifier identifier;
  unsigned numCells{0};
  std::vector<DualMemoryContainer> memoryContainer;
  std::vector<MemoryInfo> memoryInfo;

  std::unordered_map<std::type_index, std::size_t> typemap;

#ifdef ACL_DEVICE
  std::unordered_map<GraphKey, device::DeviceGraphHandle, GraphKeyHash> m_computeGraphHandles{};
  ConditionalPointersToRealsTable m_conditionalPointersToRealsTable{};
  DrConditionalPointersToRealsTable m_drConditionalPointersToRealsTable{};
  ConditionalMaterialTable m_conditionalMaterialTable{};
  ConditionalIndicesTable m_conditionalIndicesTable;
#endif

  public:
  Layer() = default;
  ~Layer() override = default;

  class CellRef {
public:
    CellRef(std::size_t id, Layer& layer, AllocationPlace place = AllocationPlace::Host)
        : layer(layer), id(id), pointers(layer.memoryInfo.size()) {
      for (std::size_t i = 0; i < layer.memoryInfo.size(); ++i) {
        if (layer.memoryInfo[i].type == MemoryType::Variable) {
          pointers[i] =
              reinterpret_cast<void*>(reinterpret_cast<char*>(layer.memoryContainer[i].get(place)) +
                                      layer.memoryInfo[i].bytes * id);
        }
      }
    }

    template <typename T>
    T& get(const Variable<T>& handle) {
      return reinterpret_cast<T*>(pointers[handle.index]);
    }

    template <typename T>
    const T& get(const Variable<T>& handle) const {
      return *reinterpret_cast<const T*>(pointers[handle.index]);
    }

    template <typename StorageT>
    typename StorageT::Type& get() {
      return *reinterpret_cast<typename StorageT::Type*>(
          pointers[layer.typemap.at(std::type_index(typeid(StorageT)))]);
    }

    template <typename StorageT>
    const typename StorageT::Type& get() const {
      return *reinterpret_cast<const typename StorageT::Type*>(
          pointers[layer.typemap.at(std::type_index(typeid(StorageT)))]);
    }

    template <typename T>
    void setPointer(const Variable<T>& handle, T* value) {
      pointers[handle.index] = reinterpret_cast<void*>(value);
    }

    template <typename T>
    T* getPointer(const Variable<T>& handle) {
      return reinterpret_cast<T*>(pointers[handle.index]);
    }

    template <typename StorageT>
    void setPointer(typename StorageT::Type* value) {
      pointers[layer.typemap.at(std::type_index(typeid(StorageT)))] =
          reinterpret_cast<void*>(value);
    }

    template <typename StorageT>
    typename StorageT::Type* getPointer() {
      return reinterpret_cast<typename StorageT::Type*>(
          pointers[layer.typemap.at(std::type_index(typeid(StorageT)))]);
    }

private:
    Layer& layer;
    std::size_t id;
    std::vector<void*> pointers;
  };

  CellRef cellRef(std::size_t id, AllocationPlace place = AllocationPlace::Host) {
    return CellRef(id, *this, place);
  }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& container : memoryContainer) {
      container.synchronizeTo(place, stream);
    }
  }

  template <typename StorageT>
  typename StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryContainer.size() > index);
    return static_cast<typename StorageT::Type*>(memoryContainer[index].get(place));
  }

  template <typename StorageT>
  void varSynchronizeTo(AllocationPlace place, void* stream) {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryContainer.size() > index);
    memoryContainer[index].synchronizeTo(place, stream);
  }

  template <typename T>
  T* var(const Variable<T>& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(memoryContainer.size() > handle.index);
    return static_cast<T*>(memoryContainer[handle.index].get(place));
  }

  template <typename T>
  void varSynchronizeTo(const Variable<T>& handle, AllocationPlace place, void* stream) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(memoryContainer.size() > handle.index);
    memoryContainer[handle.index].synchronizeTo(place, stream);
  }

  void* bucket(const Bucket& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(memoryContainer.size() > handle.index);
    return memoryContainer[handle.index].get(place);
  }

  void bucketSynchronizeTo(const Bucket& handle, AllocationPlace place, void* stream) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(memoryContainer.size() > handle.index);
    memoryContainer[handle.index].synchronizeTo(place, stream);
  }

  void* getScratchpadMemory(const ScratchpadMemory& handle,
                            AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(memoryContainer.size() > handle.index);
    return (memoryContainer[handle.index].get(place));
  }

  /// i-th bit of layerMask shall be set if data is masked on the i-th layer
  [[nodiscard]] bool isMasked(LayerMask layerMask) const {
    return (LayerMask(identifier.halo) & layerMask).any();
  }

  void setLayerType(enum LayerType layerType) { identifier.halo = layerType; }

  [[nodiscard]] enum LayerType getLayerType() const { return identifier.halo; }

  [[nodiscard]] unsigned getNumberOfCells() const { return numCells; }

  void setNumberOfCells(unsigned numberOfCells) { numCells = numberOfCells; }

  void fixPointers(const std::vector<MemoryInfo>& info,
                   const std::unordered_map<std::type_index, std::size_t>& typemap) {
    const auto count = info.size();
    memoryContainer.resize(count);
    memoryInfo.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
      memoryInfo[i] = info[i];
      memoryInfo[i].filtered = isMasked(info[i].mask) && memoryInfo[i].type == MemoryType::Variable;
    }
    this->typemap = typemap;
  }

  void setBucketSize(const Bucket& handle, size_t size) {
    assert(m_bucketSizes.size() > handle.index);
    memoryInfo[handle.index].size = size;
  }

  void setScratchpadSize(const ScratchpadMemory& handle, size_t size) {
    assert(m_scratchpadSizes.size() > handle.index);
    memoryInfo[handle.index].size = size;
  }

  size_t getBucketSize(const Bucket& handle) {
    assert(m_bucketSizes.size() > handle.index);
    return memoryInfo[handle.index].size;
  }

  void addVariableSizes(std::vector<MemoryInfo>& info, std::vector<std::size_t>& sizes) {
    for (unsigned var = 0; var < info.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Variable) {
        sizes[var] += numCells * memoryInfo[var].bytes;
      }
    }
  }

  void addBucketSizes(std::vector<MemoryInfo>& info, std::vector<std::size_t>& sizes) {
    for (unsigned var = 0; var < info.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Bucket) {
        sizes[var] += memoryInfo[var].size;
      }
    }
  }

  void findMaxScratchpadSizes(std::vector<MemoryInfo>& info, std::vector<std::size_t>& sizes) {
    for (unsigned var = 0; var < info.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Scratchpad) {
        sizes[var] = std::max(sizes[var], memoryInfo[var].size);
      }
    }
  }

  void setMemoryRegionsForVariables(const std::vector<MemoryInfo>& vars,
                                    const std::vector<DualMemoryContainer>& memory,
                                    const std::vector<size_t>& offsets) {
    for (unsigned var = 0; var < vars.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Variable) {
        memoryContainer[var].offsetFrom(memory[var], offsets[var], numCells * vars[var].bytes);
      }
    }
  }

  void setMemoryRegionsForBuckets(const std::vector<DualMemoryContainer>& memory,
                                  const std::vector<size_t>& offsets) {
    for (unsigned bucket = 0; bucket < memory.size(); ++bucket) {
      if (!memoryInfo[bucket].filtered && memoryInfo[bucket].type == MemoryType::Bucket) {
        memoryContainer[bucket].offsetFrom(
            memory[bucket], offsets[bucket], memoryInfo[bucket].size);
      }
    }
  }

  void setMemoryRegionsForScratchpads(const std::vector<DualMemoryContainer>& memory) {
    for (size_t id = 0; id < memory.size(); ++id) {
      if (!memoryInfo[id].filtered && memoryInfo[id].type == MemoryType::Scratchpad) {
        memoryContainer[id] = memory[id];
      }
    }
  }

  void touchVariables(const std::vector<MemoryInfo>& vars) {
    for (unsigned var = 0; var < vars.size(); ++var) {

      // NOTE: we don't touch device global memory because it is in a different address space
      // we will do deep-copy from the host to a device later on
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Variable &&
          (memoryContainer[var].host != nullptr)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (unsigned cell = 0; cell < numCells; ++cell) {
          memset(static_cast<char*>(memoryContainer[var].host) + cell * vars[var].bytes,
                 0,
                 vars[var].bytes);
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

#endif // SEISSOL_SRC_MEMORY_TREE_LAYER_H_
