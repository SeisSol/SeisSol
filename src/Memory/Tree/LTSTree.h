// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_
#define SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_

#include "Backmap.h"
#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Iterator.h"
#include "Common/Literals.h"
#include "Layer.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Colormap.h"
#include "Monitoring/Unit.h"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

/**
 * How the layers of a storage share their scratchpads.
 */
enum class ScratchpadSharing {
  /// All layers use the same scratchpads; no two layers may be updated concurrently.
  Shared,
  /// Each layer has scratchpads of its own.
  PerLayer
};

/*
Assigns the given value to the target object, initializing the memory in the process.

NOTE: std::copy (or the likes) do not work here, since they do not initialize the _vptr for virtual
function calls (rather, they leave it undefined), since they do merely assign `value` to `target`.
*/

template <typename T>
void initAssign(T& target, const T& value) {
  if constexpr (std::is_trivially_copyable_v<T>) {
    // if the object is trivially copyable, we may just memcpy it (it's safe to do that in this
    // case).
    std::memcpy(&target, &value, sizeof(T));
  } else {
    // otherwise, call the class/struct initializer.
    // problem: we may have an array here... So we unwrap it.
    if constexpr (std::is_array_v<T>) {
      // unwrap array, dimension by dimension...
      // example: T[N][M] yields SubT=T[M]
      using SubT = std::remove_extent_t<T>;
      auto subExtent = std::extent_v<T>;

      // for now, init element-wise... (TODO(David): we could look for something faster here, in
      // case it should ever matter)
      for (size_t i = 0; i < subExtent; ++i) {
        initAssign<SubT>(target[i], value[i]);
      }
    } else {
      // now call new here.
      new (&target) T(value);
    }
  }
  // (these two methods cannot be combined, unless we have some way for C-style arrays, i.e. S[N]
  // for <typename S, size_t N>, to use a copy constructor as well)
}

/**
  A layered storage datastructure. That is, we have some data identified by types that is clustered
  into different chunks called "layers" for historic reasons. Each layer is associated to a color
  (cf. LTSColorMap) which identifies e.g. the halo type and the LTS cluster ID it belongs to.

  The VarmapT is either GenericVarmap or SpecificVarmap.

  NOTE: currently hard-coded to the LTSColorMap type. If someone would need it, you could generalize
  that to arbitrary color maps or orderings. (we don't need that anywhere right now (late 2025))
 */
template <typename VarmapT = GenericVarmap>
class Storage {
  private:
  std::vector<DualMemoryContainer> memoryContainer_;
  std::vector<MemoryInfo> memoryInfo_;

  /// Scratchpads of the slots beyond the first one, indexed as [slot - 1][variable].
  std::vector<std::vector<DualMemoryContainer>> scratchpadSlots_;

  seissol::memory::ManagedAllocator allocator_;
  std::string name_;

  VarmapT varmap_;

  // whether the variables that change keep snapshots (see setSnapshots())
  bool snapshots_{false};

  template <typename TraitT>
  void addInternal(std::size_t index,
                   LayerMask mask,
                   size_t alignment,
                   AllocationMode allocMode,
                   bool constant,
                   std::size_t count) {
    MemoryInfo m{};
    m.alignment = alignment;
    m.mask = mask;
    m.allocMode = allocMode;
    m.constant = constant;
    m.type = TraitT::Storage;
    m.index = index;
    m.initialized = true;
    m.count = count;

    if constexpr (std::is_same_v<typename TraitT::Type, void>) {
      m.bytes = 0;
      m.bytesLayer = [count](const LayerIdentifier& identifier) {
        return dispatchConfig(identifier.config, [&](auto config) {
          using SelfT = typename TraitT::template VariantType<decltype(config)>;
          if constexpr (!std::is_same_v<void, SelfT>) {
            return sizeof(SelfT) * count;
          }
          return 0_UZ;
        });
      };
      m.alignmentLayer = [](const LayerIdentifier& identifier) {
        return dispatchConfig(identifier.config, [&](auto config) {
          using SelfT = typename TraitT::template VariantType<decltype(config)>;
          if constexpr (!std::is_same_v<void, SelfT>) {
            return alignof(SelfT);
          }
          return alignof(std::max_align_t);
        });
      };
    } else {
      using SelfT = typename TraitT::Type;
      m.bytes = sizeof(SelfT) * count;
      m.bytesLayer = [count](const LayerIdentifier& /*identifier*/) {
        return sizeof(SelfT) * count;
      };
      m.alignmentLayer = [](const LayerIdentifier& /*identifier*/) { return alignof(SelfT); };
    }

    const auto bytesLayer = m.bytesLayer;

    m.filterLayer = [mask, bytesLayer](const LayerIdentifier& identifier) {
      return mask.test(static_cast<int>(identifier.halo)) || bytesLayer(identifier) == 0;
    };

    while (memoryInfo_.size() <= index) {
      memoryInfo_.emplace_back();
    }
    memoryInfo_[index] = m;
  }

  std::optional<LTSColorMap> map_;

  std::vector<Layer<VarmapT>> layers_;

  public:
  Storage() { memoryInfo_.resize(VarmapT::MinSize); }

  ~Storage() = default;

  [[nodiscard]] const LTSColorMap& getColorMap() const { return map_.value(); }

  /// The configurations of the layers, each once, by increasing id. Every rank has the same
  /// layers, so they are the same on every rank.
  [[nodiscard]] std::vector<ConfigId> configs() const {
    std::vector<ConfigId> configs;
    configs.reserve(getColorMap().size());
    for (std::size_t color = 0; color < getColorMap().size(); ++color) {
      configs.push_back(getColorMap().argument(color).config);
    }
    std::sort(configs.begin(), configs.end());
    configs.erase(std::unique(configs.begin(), configs.end()), configs.end());
    return configs;
  }

  void setName(const std::string& name) { this->name_ = name; }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& container : memoryContainer_) {
      container.synchronizeTo(place, stream);
    }
    for (auto& slot : scratchpadSlots_) {
      for (auto& container : slot) {
        container.synchronizeTo(place, stream);
      }
    }
  }

  void setLayerCount(const LTSColorMap& map) {
    this->map_ = map;
    const auto layerCount = map.size();
    layers_.resize(layerCount);
    for (std::size_t i = 0; i < map.size(); ++i) {
      layer(i).setIdentifier(map.argument(i));
    }
  }

  void fixate() {
    memoryContainer_.resize(memoryInfo_.size());
    std::size_t id = 0;
    for (auto& leaf : this->leaves()) {
      leaf.fixPointers(id, memoryInfo_, varmap_);
      ++id;
    }
  }

  Layer<VarmapT>& layer(std::size_t index) { return layers_.at(index); }

  [[nodiscard]] const Layer<VarmapT>& layer(std::size_t index) const { return layers_.at(index); }

  Layer<VarmapT>& layer(const LayerIdentifier& id) { return layer(map_.value().colorId(id)); }

  [[nodiscard]] const Layer<VarmapT>& layer(const LayerIdentifier& id) const {
    return layer(map_.value().colorId(id));
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

  template <typename StorageT>
  const typename StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    return static_cast<typename StorageT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT>
  void varSynchronizeTo(AllocationPlace place, void* stream) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    memoryContainer_[index].synchronizeTo(place, stream);
  }

  /// Adds the memory of the variable (or bucket) `StorageT` (cf. HostData).
  template <typename StorageT>
  void varContainers(std::vector<const DualMemoryContainer*>& containers) const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer_.size() > index);
    containers.push_back(&memoryContainer_[index]);
  }

  /// Adds the memory of all variables (not of the buckets and scratchpads) (cf. HostData).
  void variableContainers(std::vector<const DualMemoryContainer*>& containers) const {
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == MemoryType::Variable) {
        containers.push_back(&memoryContainer_[var]);
      }
    }
  }

  /**
   * Whether the variables that change keep snapshots for the host, which reads them while the
   * device goes on (see DualMemoryContainer::snapshot); before allocateVariables().
   */
  void setSnapshots(bool snapshots) { snapshots_ = snapshots; }

  [[nodiscard]] bool snapshots() const { return snapshots_; }

  template <typename StorageT>
  [[nodiscard]] const MemoryInfo& info() const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryInfo_.size() > index);
    return memoryInfo_[index];
  }

  [[nodiscard]] const MemoryInfo& info(std::size_t index) const {
    assert(memoryInfo_.size() > index);
    return memoryInfo_[index];
  }

  /// The cell at `position`, seen as a cell of the configuration `Cfg`.
  template <typename Cfg>
  auto lookupRef(const StoragePosition& position, AllocationPlace place = AllocationPlace::Host) {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template cellRef<Cfg>(position.cell, place);
  }

  template <typename StorageT>
  auto& lookup(const StoragePosition& position, AllocationPlace place = AllocationPlace::Host) {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template var<StorageT>(place)[position.cell];
  }

  template <typename StorageT>
  const auto& lookup(const StoragePosition& position,
                     AllocationPlace place = AllocationPlace::Host) const {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template var<StorageT>(place)[position.cell];
  }

  /// The value of a variable at `position`, as the configuration `Cfg` holds it.
  template <typename StorageT, typename Cfg>
  auto& lookup(const Cfg& config,
               const StoragePosition& position,
               AllocationPlace place = AllocationPlace::Host) {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template var<StorageT>(config, place)[position.cell];
  }

  template <typename StorageT, typename Cfg>
  const auto& lookup(const Cfg& config,
                     const StoragePosition& position,
                     AllocationPlace place = AllocationPlace::Host) const {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template var<StorageT>(config, place)[position.cell];
  }

  [[nodiscard]] std::size_t getNumberOfVariables() const { return memoryInfo_.size(); }

  template <typename StorageT>
  void add(LayerMask mask,
           size_t alignment,
           AllocationMode allocMode,
           bool constant = false,
           std::size_t count = 1) {
    const auto index = varmap_.template add<StorageT>();
    addInternal<StorageT>(index, mask, alignment, allocMode, constant, count);
  }

  void allocateVariables() {
    std::vector<std::size_t> sizes(memoryInfo_.size());
    for (auto& leaf : this->leaves()) {
      leaf.addVariableSizes(sizes);
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == MemoryType::Variable) {
        memoryInfo_[var].size = sizes[var];
        totalSize += sizes[var];
      }
    }
    if (!name_.empty()) {
      logInfo() << "Storage" << name_ << "; variables:" << UnitByte.formatPrefix(totalSize).c_str();
    }

    std::size_t snapshotSize = 0;
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == MemoryType::Variable) {
        memoryContainer_[var].allocate(allocator_,
                                       memoryInfo_[var].size,
                                       memoryInfo_[var].alignment,
                                       memoryInfo_[var].allocMode);
        memoryContainer_[var].constant = memoryInfo_[var].constant;
        if (snapshots_) {
          memoryContainer_[var].allocateSnapshot(allocator_, memoryInfo_[var].alignment);
          if (memoryContainer_[var].snapshot != nullptr) {
            snapshotSize += memoryInfo_[var].size;
          }
        }
      }
    }
    if (!name_.empty() && snapshotSize > 0) {
      logInfo() << "Storage" << name_
                << "; snapshots of the variables:" << UnitByte.formatPrefix(snapshotSize).c_str();
    }

    std::fill(sizes.begin(), sizes.end(), 0);
    for (auto& leaf : this->leaves()) {
      leaf.setMemoryRegionsForVariables(memoryContainer_, sizes);
      leaf.addVariableSizes(sizes);
    }
  }

  void allocateBuckets() {
    std::vector<std::size_t> sizes(memoryInfo_.size());
    for (auto& leaf : this->leaves()) {
      leaf.addBucketSizes(sizes);
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == MemoryType::Bucket) {
        memoryInfo_[var].size = sizes[var];
        totalSize += sizes[var];
      }
    }
    if (!name_.empty()) {
      logInfo() << "Storage" << name_ << "; buckets:" << UnitByte.formatPrefix(totalSize).c_str();
    }

    for (std::size_t bucket = 0; bucket < memoryInfo_.size(); ++bucket) {
      if (memoryInfo_[bucket].type == MemoryType::Bucket) {
        memoryContainer_[bucket].allocate(allocator_,
                                          memoryInfo_[bucket].size,
                                          memoryInfo_[bucket].alignment,
                                          memoryInfo_[bucket].allocMode);
      }
    }

    std::fill(sizes.begin(), sizes.end(), 0);
    for (auto& leaf : this->leaves()) {
      leaf.setMemoryRegionsForBuckets(memoryContainer_, sizes);
      leaf.addBucketSizes(sizes);
    }
    touchDevice(MemoryType::Bucket);
  }

  // Allocates the scratchpads and distributes them to the leaves.
  //
  // The leaves are grouped into slots: the leaves of a slot share their scratchpads, which are
  // therefore as large as the largest demand in the slot. Leaves of the same slot must not be
  // updated concurrently; leaves of different slots may.
  void allocateScratchPads(ScratchpadSharing sharing = ScratchpadSharing::Shared) {
    const auto varCount = memoryInfo_.size();

    std::vector<std::size_t> slotOfLeaf;
    std::vector<std::vector<std::size_t>> leafSizes;
    std::size_t sizeShared = 0;
    std::size_t sizePerLayer = 0;
    {
      std::vector<std::size_t> sizesShared(varCount);
      for (auto& leaf : this->leaves()) {
        auto& sizes = leafSizes.emplace_back(varCount);
        leaf.findMaxScratchpadSizes(sizes);
        leaf.findMaxScratchpadSizes(sizesShared);
        slotOfLeaf.push_back(sharing == ScratchpadSharing::Shared ? 0 : slotOfLeaf.size());
        for (std::size_t var = 0; var < varCount; ++var) {
          if (memoryInfo_[var].type == MemoryType::Scratchpad) {
            sizePerLayer += sizes[var];
          }
        }
      }
      for (std::size_t var = 0; var < varCount; ++var) {
        if (memoryInfo_[var].type == MemoryType::Scratchpad) {
          sizeShared += sizesShared[var];
        }
      }
    }

    const auto slotCount = slotOfLeaf.empty()
                               ? std::size_t{1}
                               : *std::max_element(slotOfLeaf.begin(), slotOfLeaf.end()) + 1;

    std::vector<std::vector<std::size_t>> slotSizes(slotCount, std::vector<std::size_t>(varCount));
    for (std::size_t leaf = 0; leaf < leafSizes.size(); ++leaf) {
      auto& sizes = slotSizes[slotOfLeaf[leaf]];
      for (std::size_t var = 0; var < varCount; ++var) {
        sizes[var] = std::max(sizes[var], leafSizes[leaf][var]);
      }
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < varCount; ++var) {
      if (memoryInfo_[var].type == MemoryType::Scratchpad) {
        memoryInfo_[var].size = slotSizes[0][var];
        for (const auto& sizes : slotSizes) {
          totalSize += sizes[var];
        }
      }
    }
    if (!name_.empty()) {
      logInfo() << "Storage" << name_
                << "; scratchpads:" << UnitByte.formatPrefix(totalSize).c_str()
                << "(shared by all layers:" << UnitByte.formatPrefix(sizeShared).c_str()
                << "; one per layer:" << UnitByte.formatPrefix(sizePerLayer).c_str() << ")";
    }

    scratchpadSlots_.assign(slotCount - 1, std::vector<DualMemoryContainer>(varCount));
    const auto slotContainers = [&](std::size_t slot) -> std::vector<DualMemoryContainer>& {
      return slot == 0 ? memoryContainer_ : scratchpadSlots_[slot - 1];
    };

    for (std::size_t slot = 0; slot < slotCount; ++slot) {
      auto& containers = slotContainers(slot);
      for (std::size_t var = 0; var < varCount; ++var) {
        if (memoryInfo_[var].type == MemoryType::Scratchpad) {
          containers[var].allocate(allocator_,
                                   slotSizes[slot][var],
                                   memoryInfo_[var].alignment,
                                   memoryInfo_[var].allocMode);
        }
      }
    }

    std::size_t leafIndex = 0;
    for (auto& leaf : this->leaves()) {
      leaf.setMemoryRegionsForScratchpads(slotContainers(slotOfLeaf[leafIndex]));
      ++leafIndex;
    }
    touchDevice(MemoryType::Scratchpad);
  }

  void touchVariables() {

#pragma omp parallel
    {
      for (auto& leaf : this->leaves()) {
        leaf.touchVariables();
      }
    }
    touchDevice(MemoryType::Variable);
  }

  /**
   * Writes zeros into the device parts of the variables, buckets or scratchpads by copies (see
   * DualMemoryContainer::touchDevice()); the variables with touchVariables(), the buckets and
   * scratchpads once they are allocated.
   */
  void touchDevice(MemoryType type) {
#ifdef ACL_DEVICE
    std::vector<const DualMemoryContainer*> containers;
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == type) {
        containers.push_back(&memoryContainer_[var]);
        for (const auto& slot : scratchpadSlots_) {
          containers.push_back(&slot[var]);
        }
      }
    }
    // nothing to write (and no device to ask for the zeros) without device memory
    std::size_t largest = 0;
    for (const auto* container : containers) {
      if (container->ownDevicePart()) {
        largest = std::max(largest, container->allocationSize);
      }
    }
    if (largest == 0) {
      return;
    }
    constexpr std::size_t MaxZeros = 64UL * 1024 * 1024;
    const auto zerosSize = std::min(largest, MaxZeros);
    auto& api = device::DeviceInstance::instance().api();
    void* zeros = api.allocPinnedMem(zerosSize);
    std::memset(zeros, 0, zerosSize);
    for (const auto* container : containers) {
      container->touchDevice(zeros, zerosSize, api.getDefaultStream());
    }
    api.syncDefaultStreamWithHost();
    api.freePinnedMem(zeros);
#endif
  }

  [[nodiscard]] size_t getMaxClusterSize(LayerMask mask = LayerMask()) const {
    size_t maxClusterSize{0};
    for (const auto& leaf : this->leaves(mask)) {
      const auto currClusterSize = static_cast<size_t>(leaf.size());
      maxClusterSize = std::max(currClusterSize, maxClusterSize);
    }
    return maxClusterSize;
  }

  // helper class for the `leaves` function below; so that we can use the function like `for (auto x
  // : leaves())`
  class IteratorWrapper {
private:
    Storage<VarmapT>& node_;
    std::function<bool(const Layer<VarmapT>&)> filter_;

public:
    IteratorWrapper(Storage<VarmapT>& node, std::function<bool(const Layer<VarmapT>&)> filter)
        : node_(node), filter_(std::move(filter)) {}

    auto begin() {
      return common::FilteredIterator(node_.layers_.begin(), node_.layers_.end(), filter_);
    }

    auto end() {
      return common::FilteredIterator(node_.layers_.end(), node_.layers_.end(), filter_);
    }
  };

  // helper class for the `leaves` (const) function below; so that we can use the function like `for
  // (auto x : leaves())`
  class IteratorWrapperConst {
private:
    const Storage<VarmapT>& node_;
    std::function<bool(const Layer<VarmapT>&)> filter_;

public:
    IteratorWrapperConst(const Storage<VarmapT>& node,
                         std::function<bool(const Layer<VarmapT>&)> filter)
        : node_(node), filter_(std::move(filter)) {}

    [[nodiscard]] auto begin() const {
      return common::FilteredIterator(node_.layers_.begin(), node_.layers_.end(), filter_);
    }

    [[nodiscard]] auto end() const {
      return common::FilteredIterator(node_.layers_.end(), node_.layers_.end(), filter_);
    }
  };

  [[nodiscard]] std::size_t size(LayerMask layerMask = LayerMask()) const {
    std::size_t numCells = 0;
    for (const auto& leaf : leaves(layerMask)) {
      numCells += leaf.size();
    }
    return numCells;
  }

  template <typename F>
  IteratorWrapper filter(F&& filter) {
    return IteratorWrapper(*this, std::forward<F>(filter));
  }

  template <typename F>
  [[nodiscard]] IteratorWrapperConst filter(F&& filter) const {
    return IteratorWrapperConst(*this, std::forward<F>(filter));
  }

  IteratorWrapper leaves(LayerMask mask = LayerMask()) {
    return filter([mask](const Layer<VarmapT>& layer) {
      return !mask.test(static_cast<int>(layer.getIdentifier().halo));
    });
  }

  [[nodiscard]] IteratorWrapperConst leaves(LayerMask mask = LayerMask()) const {
    return filter([mask](const Layer<VarmapT>& layer) {
      return !mask.test(static_cast<int>(layer.getIdentifier().halo));
    });
  }

  [[nodiscard]] std::size_t numChildren() const { return layers_.size(); }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_
