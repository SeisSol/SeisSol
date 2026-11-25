// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_CORE_MEMORY_TREE_LAYER_H_
#define SEISSOL_SRC_CORE_MEMORY_TREE_LAYER_H_

#include "Core/Common/TypeId.h"
#include "Core/Memory/Batch/ConditionalTable.h"
#include "Core/Memory/Batch/DeviceGraph.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Tree/Colormap.h"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
#include <typeinfo>
#include <utils/logger.h>
#include <vector>

// TODO: remove the following

// the bits of a layer mask, used unqualified throughout
// NOLINTNEXTLINE(cppcoreguidelines-use-enum-class)
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
  DeviceOnly,
  HostDeviceCompress,
  HostDeviceCompressPinned,
  DeviceOnlyCompress,
};

enum class AllocationPlace { Host, Device };

struct DualMemoryContainer {
  void* host{nullptr};
  void* device{nullptr};
  /**
   * A copy of the device part as of the last synchronization point (see takeSnapshot()), for the
   * host to read while the device goes on, if the container keeps one (see allocateSnapshot()):
   * with separate parts, memory on the device, which gets copied to the host part later (see
   * snapshotToHost()); with unified memory, the host part itself, which is then apart from the
   * device part.
   */
  void* snapshot{nullptr};
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
    if (mode == AllocationMode::HostDeviceCompress) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::Standard);
      device = allocator.allocateMemory(
          size, alignment, seissol::memory::Memkind::DeviceGlobalCompressed);
    }
    if (mode == AllocationMode::HostDeviceCompressPinned) {
      host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::PinnedMemory);
      device = allocator.allocateMemory(
          size, alignment, seissol::memory::Memkind::DeviceGlobalCompressed);
    }
    if (mode == AllocationMode::DeviceOnlyCompress) {
      host = nullptr;
      device = allocator.allocateMemory(
          size, alignment, seissol::memory::Memkind::DeviceGlobalCompressed);
    }
    allocationMode = mode;
    allocationSize = size;
  }

  /// Whether host and device keep separate memory, which synchronizeTo() copies between.
  [[nodiscard]] bool separateParts() const {
    return allocationMode == AllocationMode::HostDeviceSplit ||
           allocationMode == AllocationMode::HostDeviceSplitPinned ||
           allocationMode == AllocationMode::HostDeviceCompress ||
           allocationMode == AllocationMode::HostDeviceCompressPinned || unifiedSnapshot();
  }

  /// Whether synchronizeTo(AllocationPlace::Host) copies the device part to the host part.
  [[nodiscard]] bool copiesToHost() const { return separateParts() && !constant; }

  /**
   * Allocates a snapshot (see `snapshot`) for data that changes, where the device can copy to the
   * host part without blocking: with separate parts, if the host part is pinned; with unified
   * memory, as a host part of its own, which the host then reads instead of the device part.
   */
  void allocateSnapshot(seissol::memory::ManagedAllocator& allocator, std::size_t alignment) {
    if (constant || allocationSize == 0 || snapshot != nullptr) {
      return;
    }
    if (allocationMode == AllocationMode::HostDeviceSplitPinned ||
        allocationMode == AllocationMode::HostDeviceCompressPinned) {
      snapshot = allocator.allocateMemory(
          allocationSize, alignment, seissol::memory::Memkind::DeviceGlobalMemory);
    }
    if (allocationMode == AllocationMode::HostDeviceUnified) {
      host = allocator.allocateMemory(
          allocationSize, alignment, seissol::memory::Memkind::DeviceUnifiedMemory);
      snapshot = host;
    }
  }

  /// Copies the device part to the snapshot, on the stream; nothing without one.
  void takeSnapshot([[maybe_unused]] void* stream) const {
#ifdef ACL_DEVICE
    if (snapshot != nullptr) {
      device::DeviceInstance::instance().api().copyBetweenAsync(
          snapshot, device, allocationSize, stream);
    }
#endif
  }

  /// Whether the snapshot still has to be copied to the host part (see snapshotToHost()).
  [[nodiscard]] bool snapshotApartFromHost() const {
    return snapshot != nullptr && snapshot != host;
  }

  /// Copies the snapshot to the host part, on the stream, where they are apart.
  void snapshotToHost([[maybe_unused]] void* stream) const {
#ifdef ACL_DEVICE
    if (snapshotApartFromHost()) {
      device::DeviceInstance::instance().api().copyFromAsync(
          host, snapshot, allocationSize, stream);
    }
#endif
  }

  void synchronizeTo(AllocationPlace place, void* stream) const {
#ifdef ACL_DEVICE
    if (unifiedSnapshot()) {
      // both parts are unified memory (see allocateSnapshot())
      if (place == AllocationPlace::Host) {
        if (!constant) {
          device::DeviceInstance::instance().api().copyBetweenAsync(
              host, device, allocationSize, stream);
        }
      } else {
        device::DeviceInstance::instance().api().copyBetweenAsync(
            device, host, allocationSize, stream);
      }
      return;
    }
    if (separateParts()) {
      if (place == AllocationPlace::Host) {
        // do not copy back constant data (we ignore the other direction for now)
        if (!constant) {
          device::DeviceInstance::instance().api().copyFromAsync(
              host, device, allocationSize, stream);
        }
      } else {
        device::DeviceInstance::instance().api().copyToAsync(device, host, allocationSize, stream);
      }
    }
    if (allocationMode == AllocationMode::HostDeviceUnified) {
      // currently broken (?)
      if (place == AllocationPlace::Host) {
        // device::DeviceInstance::instance().api().prefetchUnifiedMemTo(device::Destination::Host,
        // host, allocationSize, stream);
      } else {
        // device::DeviceInstance::instance().api().prefetchUnifiedMemTo(device::Destination::CurrentDevice,
        // device, allocationSize, stream);
      }
    }
#endif
  }

  /// Whether the device part is device memory of its own, which touchDevice() writes.
  [[nodiscard]] bool ownDevicePart() const {
    const bool deviceMemory = allocationMode == AllocationMode::HostDeviceSplit ||
                              allocationMode == AllocationMode::HostDeviceSplitPinned ||
                              allocationMode == AllocationMode::HostDeviceCompress ||
                              allocationMode == AllocationMode::HostDeviceCompressPinned ||
                              allocationMode == AllocationMode::DeviceOnly ||
                              allocationMode == AllocationMode::DeviceOnlyCompress;
    return deviceMemory && device != nullptr;
  }

  /**
   * Writes zeros into the device part by copies, as touchVariables() writes the host part, so that
   * the copy engines write it first (as compressible memory needs it on some devices, e.g. sm_80
   * and sm_86). `zeros` is a host buffer of `zerosSize` zero bytes that the device can read.
   */
  void touchDevice(const void* zeros, std::size_t zerosSize, void* stream) const {
#ifdef ACL_DEVICE
    if (!ownDevicePart()) {
      return;
    }
    for (std::size_t offset = 0; offset < allocationSize; offset += zerosSize) {
      device::DeviceInstance::instance().api().copyToAsync(
          static_cast<char*>(device) + offset,
          zeros,
          std::min(zerosSize, allocationSize - offset),
          stream);
    }
#endif
  }

  void offsetFrom(const DualMemoryContainer& container, std::size_t offset, std::size_t size) {
    host = nullptr;
    device = nullptr;
    snapshot = nullptr;

    allocationMode = container.allocationMode;
    allocationSize = size;
    if (container.host != nullptr) {
      host = reinterpret_cast<char*>(container.host) + offset;
    }
    if (container.device != nullptr) {
      device = reinterpret_cast<char*>(container.device) + offset;
    }
    if (container.snapshot != nullptr) {
      snapshot = reinterpret_cast<char*>(container.snapshot) + offset;
    }
  }

  private:
  /// Unified memory with a snapshot as a host part of its own (see allocateSnapshot()).
  [[nodiscard]] bool unifiedSnapshot() const {
    return allocationMode == AllocationMode::HostDeviceUnified && snapshot != nullptr;
  }
};

enum class MemoryType { Variable, Bucket, Scratchpad };

struct VariableDescriptor {
  static constexpr MemoryType Storage = MemoryType::Variable;
};

struct BucketDescriptor {
  static constexpr MemoryType Storage = MemoryType::Bucket;
};

struct ScratchpadDescriptor {
  static constexpr MemoryType Storage = MemoryType::Scratchpad;
};

/**
  The type a variable holds per cell in a configuration: the one it names for every configuration,
  or, where it names none (`void`), the one it names for that configuration.
 */
template <typename StorageT, typename Cfg, typename = void>
struct StorageTypeOf {
  using Type = StorageT::Type;
};

template <typename StorageT, typename Cfg>
struct StorageTypeOf<StorageT,
                     Cfg,
                     std::enable_if_t<std::is_same_v<typename StorageT::Type, void>>> {
  using Type = StorageT::template VariantType<Cfg>;
};

template <typename StorageT, typename Cfg>
using StorageType = StorageTypeOf<StorageT, Cfg>::Type;

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

template <template <typename> typename TT>
struct VariantBucket : public BucketDescriptor {
  using Type = void;

  template <typename T>
  using VariantType = TT<T>;
};

template <typename T>
struct Scratchpad : public ScratchpadDescriptor {
  using Type = T;
};

template <template <typename> typename TT>
struct VariantScratchpad : public ScratchpadDescriptor {
  using Type = void;

  template <typename T>
  using VariantType = TT<T>;
};

using FilterFunction = std::function<bool(const LayerIdentifier&)>;
using SizeFunction = std::function<std::size_t(const LayerIdentifier&)>;
using AlignmentFunction = std::function<std::size_t(const LayerIdentifier&)>;

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
  AlignmentFunction alignmentLayer;
};

template <HaloType... FilteredTypes>
bool layerFilter(const LayerIdentifier& filter) {
  return ((filter.halo == FilteredTypes) || ...);
}

/**
  Where a storage keeps its variables, buckets and scratchpads: the slot of each, an index into the
  tables of the storage and of its layers, by the id of its type (see TypeId). A storage gets its
  variables at run time (see Storage::add), so that every part of the program can add its own.
 */
class SlotMap {
  public:
  /// The slot of `StorageT`, a new one if it has none yet.
  template <typename StorageT>
  std::size_t add() {
    const auto id = TypeId<TypeDomain::Storage, StorageT>;
    if (id >= slots_.size()) {
      slots_.resize(id + 1, NoSlot);
    }
    if (slots_[id] == NoSlot) {
      slots_[id] = static_cast<std::uint32_t>(ids_.size());
      ids_.push_back(id);
    }
    return slots_[id];
  }

  /// Whether `StorageT` has a slot.
  template <typename StorageT>
  [[nodiscard]] bool has() const {
    const auto id = TypeId<TypeDomain::Storage, StorageT>;
    return id < slots_.size() && slots_[id] != NoSlot;
  }

  /// The slot of `StorageT`, which has to have one.
  template <typename StorageT>
  [[nodiscard]] std::size_t slot() const {
    if (!has<StorageT>()) [[unlikely]] {
      missing(typeid(StorageT));
    }
    return slots_[TypeId<TypeDomain::Storage, StorageT>];
  }

  /// The id of the type in the slot `slot` (see TypeId).
  [[nodiscard]] std::uint32_t idOf(std::size_t slot) const { return ids_.at(slot); }

  /// The number of slots.
  [[nodiscard]] std::size_t size() const { return ids_.size(); }

  private:
  [[noreturn]] static void missing(const std::type_info& type) {
    logError() << "The storage has no slot for" << type.name()
               << "; it has to be added (see Storage::add) before it is used.";
  }

  static constexpr std::uint32_t NoSlot = std::numeric_limits<std::uint32_t>::max();

  // by the id of the type
  std::vector<std::uint32_t> slots_;
  // by the slot
  std::vector<std::uint32_t> ids_;
};

/**
  The variables of a layer as one place (the host or the device) sees them (see Layer::view()):
  where the first cell of each is, and how many bytes a cell takes, by the id of its type (see
  TypeId). What the layer does not hold - filtered variables, buckets, scratchpads, the variables of
  other storages - is a null pointer for every cell.

  A copy can put single variables elsewhere (see relocate()), e.g. those of a cell gathered from the
  device.
 */
class LayerView {
  public:
  /// The variable `StorageT` of the cell `cell`.
  template <typename StorageT>
  [[nodiscard]] void* address(std::size_t cell) const {
    const auto id = TypeId<TypeDomain::Storage, StorageT>;
    assert(id < variables_.size());
    const auto& variable = variables_[id];
    return variable.base + variable.bytes * cell;
  }

  /// Puts the variable `StorageT` of every cell at `address`.
  template <typename StorageT>
  void relocate(void* address) {
    const auto id = TypeId<TypeDomain::Storage, StorageT>;
    assert(id < variables_.size());
    variables_[id] = Variable{.base = static_cast<char*>(address), .bytes = 0};
  }

  private:
  friend class Layer;

  struct Variable {
    char* base{nullptr};
    std::size_t bytes{0};
  };

  // by the id of the type: a table over the variables of all storages, so that a variable is found
  // with a single lookup
  std::vector<Variable> variables_;
};

/**
  The variables of a cell, for passing fewer arguments (equivalent to the `LocalData` /
  `NeighborData` in older SeisSol). It holds no more than the variables of its layer as one place
  sees them (see LayerView) and the cell: a variable is found when it is asked for, so that making a
  reference costs nothing.

  The reference is to a cell of a configuration `Cfg`: a variable that holds a type per
  configuration is seen as the one of `Cfg`.

  Example: if you have a registered type Tx, then the following are equivalent:

  const auto& data = layer.var<Tx>(Cfg(), AllocationPlace::Device)[123];

  to:

  const auto ref = layer.cellRef<Cfg>(123, AllocationPlace::Device);
  const auto& data = ref.get<Tx>();
 */
template <typename Cfg>
class CellRef {
  public:
  CellRef(const LayerView& view, std::size_t cell) : view_(&view), cell_(cell) {}

  template <typename StorageT>
  StorageType<StorageT, Cfg>& get() {
    return *static_cast<StorageType<StorageT, Cfg>*>(view_->template address<StorageT>(cell_));
  }

  template <typename StorageT>
  [[nodiscard]] const StorageType<StorageT, Cfg>& get() const {
    return *static_cast<const StorageType<StorageT, Cfg>*>(
        view_->template address<StorageT>(cell_));
  }

  /// The view of the layer of the cell, e.g. to derive another one (see LayerView::relocate()).
  [[nodiscard]] const LayerView& view() const { return *view_; }

  [[nodiscard]] std::size_t cell() const { return cell_; }

  private:
  const LayerView* view_;
  std::size_t cell_;
};

/**
  A layer (i.e. contiguous array of data) in a Storage structure.
 */
class Layer {
  private:
  LayerIdentifier identifier_;
  std::size_t numCells_{0};
  std::vector<DualMemoryContainer> memoryContainer_;
  std::vector<MemoryInfo> memoryInfo_;

  std::size_t posId_{0};

  SlotMap slots_;

  // the variables as the host and as the device see them (see view())
  std::array<LayerView, 2> views_;

#ifdef ACL_DEVICE
  std::unordered_map<GraphKey, device::DeviceGraphHandle, GraphKeyHash> computeGraphHandles_;
#endif

  recording::ConditionalTable conditionalTable_;

  static std::size_t placeIndex(AllocationPlace place) {
    return place == AllocationPlace::Host ? 0 : 1;
  }

  /// Finds the variables anew for the views, once their memory is known.
  void updateViews() {
    const auto ids = typeCount(TypeDomain::Storage);
    for (const auto place : {AllocationPlace::Host, AllocationPlace::Device}) {
      auto& variables = views_[placeIndex(place)].variables_;
      variables.assign(ids, LayerView::Variable{});
      for (std::size_t slot = 0; slot < memoryInfo_.size(); ++slot) {
        const auto& info = memoryInfo_[slot];
        auto* base = static_cast<char*>(memoryContainer_[slot].get(place));
        if (info.type == MemoryType::Variable && info.initialized && !info.filtered &&
            base != nullptr) {
          variables[slots_.idOf(slot)] = LayerView::Variable{.base = base, .bytes = info.bytes};
        }
      }
    }
  }

  public:
  Layer() = default;
  ~Layer() = default;

  [[nodiscard]] std::size_t id() const { return posId_; }

  /// The variables of the layer as `place` sees them.
  [[nodiscard]] const LayerView& view(AllocationPlace place = AllocationPlace::Host) const {
    return views_[placeIndex(place)];
  }

  /// The cell `id`, as a cell of the configuration `Cfg` (see CellRef).
  template <typename Cfg>
  CellRef<Cfg> cellRef(std::size_t id, AllocationPlace place = AllocationPlace::Host) {
    return CellRef<Cfg>(view(place), id);
  }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (const auto& container : memoryContainer_) {
      container.synchronizeTo(place, stream);
    }
  }

  void* varUntyped(std::size_t index, AllocationPlace place = AllocationPlace::Host) {
    assert(index != std::numeric_limits<std::size_t>::max());
    assert(memoryContainer_.size() > index);
    return memoryContainer_[index].get(place);
  }

  template <typename StorageT>
  StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) {
    const auto index = slots_.slot<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<StorageT::Type*>(memoryContainer_[index].get(place));
  }

  /// The values of a variable as the configuration `ConfigT` holds them.
  template <typename StorageT, typename ConfigT>
  StorageType<StorageT, ConfigT>* var(const ConfigT& /*config*/,
                                      AllocationPlace place = AllocationPlace::Host)
    requires(!std::is_same_v<ConfigT, AllocationPlace>)
  {
    const auto index = slots_.slot<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<StorageType<StorageT, ConfigT>*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT>
  [[nodiscard]] const StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) const {
    const auto index = slots_.slot<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<StorageT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT, typename ConfigT>
  [[nodiscard]] const StorageType<StorageT, ConfigT>*
      var(const ConfigT& /*config*/, AllocationPlace place = AllocationPlace::Host) const
    requires(!std::is_same_v<ConfigT, AllocationPlace>)
  {
    const auto index = slots_.slot<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<StorageType<StorageT, ConfigT>*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT>
  void varSynchronizeTo(AllocationPlace place, void* stream) {
    const auto index = slots_.slot<StorageT>();
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

  void fixPointers(std::size_t id, const std::vector<MemoryInfo>& info, const SlotMap& slots) {
    posId_ = id;
    const auto count = info.size();
    memoryContainer_.resize(count);
    memoryInfo_.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
      memoryInfo_[i] = info[i];
      if (memoryInfo_[i].initialized) {
        memoryInfo_[i].filtered = info[i].filterLayer(identifier_);
        memoryInfo_[i].bytes = info[i].bytesLayer(identifier_);

        // alignment sanity check
        const auto neededAlignment = memoryInfo_[i].alignmentLayer(identifier_);
        if (neededAlignment > memoryInfo_[i].alignment) {
          logError() << "Alignment mismatch:" << memoryInfo_[i].alignment << "given; needed"
                     << neededAlignment;
        }
      }
    }
    slots_ = slots;
    updateViews();
  }

  template <typename StorageT>
  void setEntrySize(size_t size) {
    const auto index = slots_.slot<StorageT>();
    assert(memoryInfo_.size() > index);
    static_assert(StorageT::Storage == MemoryType::Bucket ||
                  StorageT::Storage == MemoryType::Scratchpad);
    memoryInfo_[index].size = size;
  }

  template <typename StorageT>
  size_t getEntrySize() {
    const auto index = slots_.slot<StorageT>();
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
    updateViews();
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

#pragma omp for schedule(static) nowait
        for (std::size_t cell = 0; cell < numCells_; ++cell) {
          auto* cellPointer =
              static_cast<char*>(memoryContainer_[var].host) + cell * memoryInfo_[var].bytes;
          memset(cellPointer, 0, memoryInfo_[var].bytes);
        }
      }
    }
  }

  /// The batches recorded for the layer, by their keys (see recording::ConditionalKey).
  recording::ConditionalTable& getConditionalTable() { return conditionalTable_; }

  [[nodiscard]] const recording::ConditionalTable& getConditionalTable() const {
    return conditionalTable_;
  }

#ifdef ACL_DEVICE
  device::DeviceGraphHandle getDeviceComputeGraphHandle(GraphKey graphKey) {
    if (computeGraphHandles_.contains(graphKey)) {
      return computeGraphHandles_[graphKey];
    } else {
      return device::DeviceGraphHandle();
    }
  }

  void updateDeviceComputeGraphHandle(GraphKey graphKey,
                                      const device::DeviceGraphHandle& graphHandle) {
    assert(!computeGraphHandles_.contains(graphKey) &&
           "an entry of hash table must be empty on write");
    if (graphHandle.isInitialized()) {
      computeGraphHandles_[graphKey] = graphHandle;
    }
  }
#endif // ACL_DEVICE
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_CORE_MEMORY_TREE_LAYER_H_
