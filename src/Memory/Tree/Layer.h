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
#include "Memory/Tree/Colormap.h"

#include <bitset>
#include <cstring>
#include <limits>
#include <type_traits>
#include <typeindex>
#include <utils/logger.h>
#include <variant>
#include <vector>

// TODO: remove the following

// NOLINTBEGIN (-cppcoreguidelines-use-enum-class)
enum LayerType { Ghost = (1U << 0U), Copy = (1U << 1U), Interior = (1U << 2U), NumLayers = 3U };
// NOLINTEND (-cppcoreguidelines-use-enum-class)

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
  void* host{nullptr};
  void* device{nullptr};
  AllocationMode allocationMode{AllocationMode::HostOnly};
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

  void synchronizeTo(AllocationPlace place, void* stream) const {
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
  std::shared_ptr<int> handle{std::make_shared<int>()};

  [[nodiscard]] int* pointer() const { return handle.get(); }

  MemoryHandle() = default;
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

using FilterFunction = std::function<bool(const LayerIdentifier&)>;
using SizeFunction = std::function<std::size_t(const LayerIdentifier&)>;

struct MemoryInfo {
  size_t index{};
  size_t bytes{};
  size_t alignment{};
  size_t size{};
  size_t count{};
  LayerMask mask;
  AllocationMode allocMode{AllocationMode::HostOnly};
  MemoryType type{MemoryType::Variable};
  bool constant{false};
  bool filtered{false};
  bool initialized{false};

  SizeFunction bytesLayer;
  FilterFunction filterLayer;
};

template <HaloType... FilteredTypes>
bool layerFilter(const LayerIdentifier& filter) {
  return ((filter.halo == FilteredTypes) || ...);
}

/**
  A generic Varmap which does not have any type restrictions.
  Slower in access than a SpecificVarmap, but it doesn't require any pre-defined type lists to work.
 */
struct GenericVarmap {
  std::unordered_map<std::type_index, std::size_t> typemap;
  std::unordered_map<int*, std::size_t> handlemap;
  std::size_t count = 0;

  template <typename T>
  std::size_t add() {
    const auto idx = count;
    typemap[std::type_index(typeid(T))] = idx;
    ++count;
    return idx;
  }

  template <typename T>
  std::size_t add(T& handle) {
    const auto idx = count;
    handlemap[handle.pointer()] = idx;
    ++count;
    return idx;
  }

  template <typename T>
  [[nodiscard]] std::size_t index() const {
    return typemap.at(std::type_index(typeid(T)));
  }

  template <typename T>
  [[nodiscard]] std::size_t index(T& handle) const {
    return handlemap.at(handle.pointer());
  }

  static constexpr std::size_t MinSize = 0;

  [[nodiscard]] std::vector<void*> pointerContainer() const { return std::vector<void*>(count); }

  using PointerContainerT = std::vector<void*>;
};

/**
  A Varmap with a type list. It is faster in accessing, but needs all types pre-registered (leading
  to usually veeeeery long type lists; try to subclass to avoid the error messages being too long).
  */
template <typename... Types>
struct SpecificVarmap {
  template <typename T, typename Thead, typename... Ttail>
  [[nodiscard]] std::size_t constexpr innerIndex(std::size_t value = 0) const {
    if constexpr (std::is_same_v<T, Thead>) {
      return value;
    } else if constexpr (sizeof...(Ttail) == 0) {
      static_assert(sizeof(T) == 0, "Type not found.");
      return std::numeric_limits<std::size_t>::max();
    } else {
      return innerIndex<T, Ttail...>(value + 1);
    }
  }

  template <typename T>
  std::size_t add() {
    return index<T>();
  }

  template <typename T>
  std::size_t add(T& handle) {
    return index(handle);
  }

  template <typename T>
  [[nodiscard]] std::size_t index() const {
    return innerIndex<T, Types...>();
  }

  template <typename T>
  [[nodiscard]] std::size_t index(T& /*handle*/) const {
    static_assert(sizeof(T) == 0, "Type not found.");
    return 0;
  }

  static constexpr std::size_t MinSize = sizeof...(Types);

  using PointerContainerT = std::array<void*, MinSize>;

  [[nodiscard]] std::array<void*, MinSize> pointerContainer() const {
    return std::array<void*, MinSize>();
  }
};

/**
  A layer (i.e. contiguous array of data) in a Storage structure.
 */
template <typename VarmapT = GenericVarmap>
class Layer {
  private:
  LayerIdentifier identifier_;
  std::size_t numCells_{0};
  std::vector<DualMemoryContainer> memoryContainer_;
  std::vector<MemoryInfo> memoryInfo_;

  std::size_t posId_{0};

  VarmapT varmap_;

#ifdef ACL_DEVICE
  std::unordered_map<GraphKey, device::DeviceGraphHandle, GraphKeyHash> computeGraphHandles_;
#endif

  recording::ConditionalPointersToRealsTable conditionalPointersToRealsTable_;
  recording::DrConditionalPointersToRealsTable drConditionalPointersToRealsTable_;
  recording::ConditionalMaterialTable conditionalMaterialTable_;
  recording::ConditionalIndicesTable conditionalIndicesTable_;

  public:
  Layer() = default;
  ~Layer() = default;

  [[nodiscard]] std::size_t id() const { return posId_; }

  /**
    Contains references to a specific cell for all stored internal data.
    Useful e.g. for passing less arguments.

    (equivalent to the `LocalData` / `NeighborData` in older SeisSol)

    Example: if you have a registered type Tx, then the following are equivalent:

    const auto& data = layer.var<Tx>(AllocationPlace::Device)[123];

    to:

    const auto ref = layer.cellRef(123, AllocationPlace::Device);
    const auto& data = ref.get<Tx>();
   */
  class CellRef {
public:
    CellRef(std::size_t id,
            const VarmapT& varmap,
            Layer& layer,
            AllocationPlace place = AllocationPlace::Host)
        : pointers_(varmap.pointerContainer()) {
      for (std::size_t i = 0; i < layer.memoryInfo_.size(); ++i) {
        if (layer.memoryInfo_[i].type == MemoryType::Variable && !layer.memoryInfo_[i].filtered &&
            layer.memoryInfo_[i].initialized) {
          pointers_[i] = reinterpret_cast<void*>(
              reinterpret_cast<char*>(layer.memoryContainer_[i].get(place)) +
              layer.memoryInfo_[i].bytes * id);
        } else {
          pointers_[i] = nullptr;
        }
      }
    }

    template <typename HandleT>
    typename HandleT::Type& get(const HandleT& handle) {
      return *reinterpret_cast<typename HandleT::Type*>(pointers_[varmap_.index(handle)]);
    }

    template <typename HandleT>
    [[nodiscard]] const typename HandleT::Type& get(const HandleT& handle) const {
      return *reinterpret_cast<const typename HandleT::Type*>(pointers_[varmap_.index(handle)]);
    }

    template <typename StorageT>
    typename StorageT::Type& get() {
      return *reinterpret_cast<typename StorageT::Type*>(
          pointers_[varmap_.template index<StorageT>()]);
    }

    template <typename StorageT>
    [[nodiscard]] const typename StorageT::Type& get() const {
      return *reinterpret_cast<const typename StorageT::Type*>(
          pointers_[varmap_.template index<StorageT>()]);
    }

    template <typename HandleT>
    void setPointer(const HandleT& handle, typename HandleT::Type* value) {
      pointers_[varmap_.index(handle)] = reinterpret_cast<void*>(value);
    }

    template <typename HandleT>
    typename HandleT::Type* getPointer(const HandleT& handle) {
      return reinterpret_cast<typename HandleT::Type*>(pointers_[varmap_.index(handle)]);
    }

    template <typename StorageT>
    void setPointer(typename StorageT::Type* value) {
      pointers_[varmap_.template index<StorageT>()] = reinterpret_cast<void*>(value);
    }

    template <typename StorageT>
    typename StorageT::Type* getPointer() {
      return reinterpret_cast<typename StorageT::Type*>(
          pointers_[varmap_.template index<StorageT>()]);
    }

private:
    VarmapT varmap_;
    typename VarmapT::PointerContainerT pointers_;
  };

  CellRef cellRef(std::size_t id, AllocationPlace place = AllocationPlace::Host) {
    return CellRef(id, varmap_, *this, place);
  }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& container : memoryContainer_) {
      container.synchronizeTo(place, stream);
    }
  }

  void* varUntyped(std::size_t index, AllocationPlace place = AllocationPlace::Host) {
    assert(index != std::numeric_limits<std::size_t>::max());
    assert(memoryContainer_.size() > index);
    return memoryContainer_[index].get(place);
  }

  template <typename StorageT>
  typename StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<typename StorageT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT, typename ConfigT>
  typename StorageT::template VariantType<ConfigT>*
      var(const ConfigT& /*...*/, AllocationPlace place = AllocationPlace::Host) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<typename StorageT::template VariantType<ConfigT>*>(
        memoryContainer_[index].get(place));
  }

  template <typename StorageT>
  [[nodiscard]] const typename StorageT::Type*
      var(AllocationPlace place = AllocationPlace::Host) const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<typename StorageT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT, typename ConfigT>
  [[nodiscard]] const typename StorageT::template VariantType<ConfigT>*
      var(const ConfigT& /*...*/, AllocationPlace place = AllocationPlace::Host) const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<typename StorageT::template VariantType<ConfigT>*>(
        memoryContainer_[index].get(place));
  }

  template <typename StorageT>
  void varSynchronizeTo(AllocationPlace place, void* stream) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    memoryContainer_[index].synchronizeTo(place, stream);
  }

  template <typename HandleT>
  typename HandleT::Type* var(const HandleT& handle,
                              AllocationPlace place = AllocationPlace::Host) {
    const auto index = varmap_.index(handle);
    assert(memoryContainer_.size() > index);
    return static_cast<typename HandleT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename HandleT, typename ConfigT>
  typename HandleT::template VariantType<ConfigT>*
      var(const HandleT& handle,
          const ConfigT& /*...*/,
          AllocationPlace place = AllocationPlace::Host) {
    const auto index = varmap_.index(handle);
    assert(memoryContainer_.size() > index);
    return static_cast<typename HandleT::template VariantType<ConfigT>*>(
        memoryContainer_[index].get(place));
  }

  template <typename HandleT>
  [[nodiscard]] const typename HandleT::Type*
      var(const HandleT& handle, AllocationPlace place = AllocationPlace::Host) const {
    const auto index = varmap_.index(handle);
    assert(memoryContainer_.size() > index);
    return static_cast<typename HandleT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename HandleT, typename ConfigT>
  [[nodiscard]] const typename HandleT::template VariantType<ConfigT>*
      var(const HandleT& handle,
          const ConfigT& /*...*/,
          AllocationPlace place = AllocationPlace::Host) const {
    const auto index = varmap_.index(handle);
    assert(memoryContainer_.size() > index);
    return static_cast<typename HandleT::template VariantType<ConfigT>*>(
        memoryContainer_[index].get(place));
  }

  template <typename HandleT>
  void varSynchronizeTo(const HandleT& handle, AllocationPlace place, void* stream) {
    const auto index = varmap_.index(handle);
    assert(memoryContainer_.size() > index);
    memoryContainer_[index].synchronizeTo(place, stream);
  }

  /// i-th bit of layerMask shall be set if data is masked on the i-th layer
  [[nodiscard]] bool isMasked(LayerMask layerMask) const {
    return layerMask.test(static_cast<int>(getIdentifier().halo));
  }

  [[nodiscard]] const LayerIdentifier& getIdentifier() const { return identifier_; }

  void setIdentifier(const LayerIdentifier& identifier) { this->identifier_ = identifier; }

  [[nodiscard]] std::size_t size() const { return numCells_; }

  void setNumberOfCells(std::size_t numberOfCells) { numCells_ = numberOfCells; }

  void fixPointers(std::size_t id, const std::vector<MemoryInfo>& info, const VarmapT& varmap) {
    posId_ = id;
    const auto count = info.size();
    memoryContainer_.resize(count);
    memoryInfo_.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
      memoryInfo_[i] = info[i];
      if (memoryInfo_[i].initialized) {
        memoryInfo_[i].filtered = info[i].filterLayer(identifier_);
        memoryInfo_[i].bytes = info[i].bytesLayer(identifier_);
      }
    }
    this->varmap_ = varmap;
  }

  template <typename HandleT>
  void setEntrySize(const HandleT& handle, size_t size) {
    const auto index = varmap_.index(handle);
    assert(memoryInfo_.size() > index);
    static_assert(HandleT::Storage == MemoryType::Bucket ||
                  HandleT::Storage == MemoryType::Scratchpad);
    memoryInfo_[index].size = size;
  }

  template <typename HandleT>
  size_t getEntrySize(const HandleT& handle) {
    const auto index = varmap_.index(handle);
    assert(memoryInfo_.size() > index);
    return memoryInfo_[index].size;
  }

  template <typename StorageT>
  void setEntrySize(size_t size) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryInfo_.size() > index);
    static_assert(StorageT::Storage == MemoryType::Bucket ||
                  StorageT::Storage == MemoryType::Scratchpad);
    memoryInfo_[index].size = size;
  }

  template <typename StorageT>
  size_t getEntrySize() {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryInfo_.size() > index);
    return memoryInfo_[index].size;
  }

  void addVariableSizes(std::vector<std::size_t>& sizes) {
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (!memoryInfo_[var].filtered && memoryInfo_[var].type == MemoryType::Variable) {
        sizes[var] += numCells_ * memoryInfo_[var].bytes;
      }
    }
  }

  void addBucketSizes(std::vector<std::size_t>& sizes) {
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (!memoryInfo_[var].filtered && memoryInfo_[var].type == MemoryType::Bucket) {
        sizes[var] += memoryInfo_[var].size;
      }
    }
  }

  void findMaxScratchpadSizes(std::vector<std::size_t>& sizes) {
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (!memoryInfo_[var].filtered && memoryInfo_[var].type == MemoryType::Scratchpad) {
        sizes[var] = std::max(sizes[var], memoryInfo_[var].size);
      }
    }
  }

  void setMemoryRegionsForVariables(const std::vector<DualMemoryContainer>& memory,
                                    const std::vector<size_t>& offsets) {
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (!memoryInfo_[var].filtered && memoryInfo_[var].type == MemoryType::Variable) {
        memoryContainer_[var].offsetFrom(
            memory[var], offsets[var], numCells_ * memoryInfo_[var].bytes);
      }
    }
  }

  void setMemoryRegionsForBuckets(const std::vector<DualMemoryContainer>& memory,
                                  const std::vector<size_t>& offsets) {
    for (std::size_t bucket = 0; bucket < memory.size(); ++bucket) {
      if (!memoryInfo_[bucket].filtered && memoryInfo_[bucket].type == MemoryType::Bucket) {
        memoryContainer_[bucket].offsetFrom(
            memory[bucket], offsets[bucket], memoryInfo_[bucket].size);
      }
    }
  }

  void setMemoryRegionsForScratchpads(const std::vector<DualMemoryContainer>& memory) {
    for (size_t id = 0; id < memory.size(); ++id) {
      if (!memoryInfo_[id].filtered && memoryInfo_[id].type == MemoryType::Scratchpad) {
        memoryContainer_[id] = memory[id];
      }
    }
  }

  void touchVariables() {
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {

      // NOTE: we don't touch device global memory because it is in a different address space
      // we will do deep-copy from the host to a device later on
      if (!memoryInfo_[var].filtered && memoryInfo_[var].type == MemoryType::Variable &&
          (memoryContainer_[var].host != nullptr) && numCells_ > 0) {
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
        for (std::size_t cell = 0; cell < numCells_; ++cell) {
          auto* cellPointer =
              static_cast<char*>(memoryContainer_[var].host) + cell * memoryInfo_[var].bytes;
          memset(cellPointer, 0, memoryInfo_[var].bytes);
        }
      }
    }
  }

  template <typename F>
  void wrap(F&& function) {
    std::visit(std::forward<F>(function), identifier_.config);
  }

  template <typename InnerKeyType>
  auto& getConditionalTable() {
    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Wp>) {
      return conditionalPointersToRealsTable_;
    }

    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Dr>) {
      return drConditionalPointersToRealsTable_;
    }

    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Material>) {
      return conditionalMaterialTable_;
    }

    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Indices>) {
      return conditionalIndicesTable_;
    }
  }

  template <typename InnerKeyType>
  const auto& getConditionalTable() const {
    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Wp>) {
      return conditionalPointersToRealsTable_;
    }

    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Dr>) {
      return drConditionalPointersToRealsTable_;
    }

    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Material>) {
      return conditionalMaterialTable_;
    }

    if constexpr (std::is_same_v<InnerKeyType, recording::inner_keys::Indices>) {
      return conditionalIndicesTable_;
    }
  }

#ifdef ACL_DEVICE
  device::DeviceGraphHandle getDeviceComputeGraphHandle(GraphKey graphKey) {
    if (computeGraphHandles_.find(graphKey) != computeGraphHandles_.end()) {
      return computeGraphHandles_[graphKey];
    } else {
      return device::DeviceGraphHandle();
    }
  }

  void updateDeviceComputeGraphHandle(GraphKey graphKey, device::DeviceGraphHandle graphHandle) {
    assert(computeGraphHandles_.find(graphKey) == computeGraphHandles_.end() &&
           "an entry of hash table must be empty on write");
    if (graphHandle.isInitialized()) {
      computeGraphHandles_[graphKey] = graphHandle;
    }
  }
#endif // ACL_DEVICE
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_TREE_LAYER_H_
