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

enum class MemoryType { Variable, Bucket, Scratchpad };

struct MemoryHandle {
  std::shared_ptr<int> handle;

  [[nodiscard]] int* pointer() const { return handle.get(); }

  MemoryHandle() : handle(std::make_shared<int>()) {}
};

struct VariableDescriptor : public MemoryHandle {
  static constexpr MemoryType Storage = MemoryType::Variable;
};

struct BucketDescriptor : public MemoryHandle {
  static constexpr MemoryType Storage = MemoryType::Bucket;
};

struct ScratchpadDescriptor : public MemoryHandle {
  static constexpr MemoryType Storage = MemoryType::Scratchpad;
};

template <typename T>
struct Variable : public VariableDescriptor {
  using Type = T;
};

template <template <typename> typename TT>
struct VariantVariable : public VariableDescriptor {
  using Type = void;

  template <typename T>
  using VariantType = TT<T>;
};

template <typename T>
struct Bucket : public BucketDescriptor {
  using Type = T;
};

template <typename T>
struct Scratchpad : public ScratchpadDescriptor {
  using Type = T;
};

struct LayerIdentifier {
  enum LayerType halo;
  int lts;
  std::variant<Config> config;
};

using FilterFunction = std::function<bool(const LayerIdentifier&)>;
using SizeFunction = std::function<std::size_t(const LayerIdentifier&)>;

struct MemoryInfo {
  size_t index{};
  size_t bytes{};
  size_t alignment{};
  size_t size{};
  LayerMask mask;
  AllocationMode allocMode;
  MemoryType type;
  bool constant{false};
  bool filtered{false};

  SizeFunction bytesLayer;
  FilterFunction filterLayer;
};

template <LayerType... FilteredTypes>
bool layerFilter(const LayerIdentifier& filter) {
  return ((filter.halo == FilteredTypes) || ...);
}

class Layer : public Node {
  private:
  LayerIdentifier identifier;
  std::size_t numCells{0};
  std::vector<DualMemoryContainer> memoryContainer;
  std::vector<MemoryInfo> memoryInfo;

  std::unordered_map<std::type_index, std::size_t> typemap;
  std::unordered_map<int*, std::size_t> handlemap;

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
        : layer(layer), pointers(layer.memoryInfo.size()) {
      for (std::size_t i = 0; i < layer.memoryInfo.size(); ++i) {
        if (layer.memoryInfo[i].type == MemoryType::Variable && !layer.memoryInfo[i].filtered) {
          pointers[i] =
              reinterpret_cast<void*>(reinterpret_cast<char*>(layer.memoryContainer[i].get(place)) +
                                      layer.memoryInfo[i].bytes * id);
        }
      }
    }

    template <typename HandleT>
    typename HandleT::Type& get(const HandleT& handle) {
      return *reinterpret_cast<typename HandleT::Type*>(
          pointers[layer.handlemap.at(handle.pointer())]);
    }

    template <typename HandleT>
    const typename HandleT::Type& get(const HandleT& handle) const {
      return *reinterpret_cast<const typename HandleT::Type*>(
          pointers[layer.handlemap.at(handle.pointer())]);
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

    template <typename HandleT>
    void setPointer(const HandleT& handle, typename HandleT::Type* value) {
      pointers[layer.handlemap.at(handle.pointer())] = reinterpret_cast<void*>(value);
    }

    template <typename HandleT>
    typename HandleT::Type* getPointer(const HandleT& handle) {
      return reinterpret_cast<typename HandleT::Type*>(
          pointers[layer.handlemap.at(handle.pointer())]);
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

  template <typename StorageT, typename ConfigT>
  typename StorageT::template VariantType<ConfigT>*
      var(const ConfigT& /*...*/, AllocationPlace place = AllocationPlace::Host) {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryContainer.size() > index);
    return static_cast<typename StorageT::template VariantType<ConfigT>*>(
        memoryContainer[index].get(place));
  }

  template <typename StorageT>
  void varSynchronizeTo(AllocationPlace place, void* stream) {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryContainer.size() > index);
    memoryContainer[index].synchronizeTo(place, stream);
  }

  template <typename HandleT>
  typename HandleT::Type* var(const HandleT& handle,
                              AllocationPlace place = AllocationPlace::Host) {
    const auto index = handlemap.at(handle.pointer());
    assert(memoryContainer.size() > index);
    return static_cast<typename HandleT::Type*>(memoryContainer[index].get(place));
  }

  template <typename HandleT, typename ConfigT>
  typename HandleT::template VariantType<ConfigT>*
      var(const HandleT& handle,
          const ConfigT& /*...*/,
          AllocationPlace place = AllocationPlace::Host) {
    const auto index = handlemap.at(handle.pointer());
    assert(memoryContainer.size() > index);
    return static_cast<typename HandleT::template VariantType<ConfigT>*>(
        memoryContainer[index].get(place));
  }

  template <typename HandleT>
  void varSynchronizeTo(const HandleT& handle, AllocationPlace place, void* stream) {
    const auto index = handlemap.at(handle.pointer());
    assert(memoryContainer.size() > index);
    memoryContainer[index].synchronizeTo(place, stream);
  }

  /// i-th bit of layerMask shall be set if data is masked on the i-th layer
  [[nodiscard]] bool isMasked(LayerMask layerMask) const {
    return (LayerMask(identifier.halo) & layerMask).any();
  }

  void setLayerType(enum LayerType layerType) { identifier.halo = layerType; }

  [[nodiscard]] const LayerIdentifier& getIdentifier() const { return identifier; }

  [[nodiscard]] std::size_t size() const { return numCells; }

  void setNumberOfCells(std::size_t numberOfCells) { numCells = numberOfCells; }

  void fixPointers(const std::vector<MemoryInfo>& info,
                   const std::unordered_map<std::type_index, std::size_t>& typemap,
                   const std::unordered_map<int*, std::size_t>& handlemap) {
    const auto count = info.size();
    memoryContainer.resize(count);
    memoryInfo.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
      memoryInfo[i] = info[i];
      memoryInfo[i].filtered = info[i].filterLayer(identifier);
      memoryInfo[i].bytes = info[i].bytesLayer(identifier);
    }
    this->typemap = typemap;
    this->handlemap = handlemap;
  }

  template <typename HandleT>
  void setEntrySize(const HandleT& handle, size_t size) {
    const auto index = handlemap.at(handle.pointer());
    assert(memoryInfo.size() > index);
    static_assert(HandleT::Storage == MemoryType::Bucket ||
                  HandleT::Storage == MemoryType::Scratchpad);
    memoryInfo[index].size = size;
  }

  template <typename HandleT>
  size_t getEntrySize(const HandleT& handle) {
    const auto index = handlemap.at(handle.pointer());
    assert(memoryInfo.size() > index);
    return memoryInfo[index].size;
  }

  template <typename StorageT>
  void setEntrySize(size_t size) {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryInfo.size() > index);
    static_assert(StorageT::Storage == MemoryType::Bucket ||
                  StorageT::Storage == MemoryType::Scratchpad);
    memoryInfo[index].size = size;
  }

  template <typename StorageT>
  size_t getEntrySize() {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryInfo.size() > index);
    return memoryInfo[index].size;
  }

  void addVariableSizes(std::vector<MemoryInfo>& info, std::vector<std::size_t>& sizes) {
    for (std::size_t var = 0; var < info.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Variable) {
        sizes[var] += numCells * memoryInfo[var].bytes;
      }
    }
  }

  void addBucketSizes(std::vector<MemoryInfo>& info, std::vector<std::size_t>& sizes) {
    for (std::size_t var = 0; var < info.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Bucket) {
        sizes[var] += memoryInfo[var].size;
      }
    }
  }

  void findMaxScratchpadSizes(std::vector<MemoryInfo>& info, std::vector<std::size_t>& sizes) {
    for (std::size_t var = 0; var < info.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Scratchpad) {
        sizes[var] = std::max(sizes[var], memoryInfo[var].size);
      }
    }
  }

  void setMemoryRegionsForVariables(const std::vector<MemoryInfo>& vars,
                                    const std::vector<DualMemoryContainer>& memory,
                                    const std::vector<size_t>& offsets) {
    for (std::size_t var = 0; var < vars.size(); ++var) {
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Variable) {
        memoryContainer[var].offsetFrom(memory[var], offsets[var], numCells * vars[var].bytes);
      }
    }
  }

  void setMemoryRegionsForBuckets(const std::vector<DualMemoryContainer>& memory,
                                  const std::vector<size_t>& offsets) {
    for (std::size_t bucket = 0; bucket < memory.size(); ++bucket) {
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
    for (std::size_t var = 0; var < vars.size(); ++var) {

      // NOTE: we don't touch device global memory because it is in a different address space
      // we will do deep-copy from the host to a device later on
      if (!memoryInfo[var].filtered && memoryInfo[var].type == MemoryType::Variable &&
          (memoryContainer[var].host != nullptr) && numCells > 0) {
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
        for (std::size_t cell = 0; cell < numCells; ++cell) {
          auto* cellPointer =
              static_cast<char*>(memoryContainer[var].host) + cell * vars[var].bytes;
          memset(cellPointer, 0, vars[var].bytes);
        }
      }
    }
  }

  template <typename F>
  void wrap(F&& function) {
    std::visit(std::forward<F>(function), identifier.config);
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
