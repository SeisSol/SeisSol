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
#include "Common/Iterator.h"
#include "Config.h"
#include "Layer.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/Colormap.h"
#include "Monitoring/Unit.h"

#include <type_traits>
#include <utility>
#include <utils/logger.h>

namespace seissol::initializer {

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

  seissol::memory::ManagedAllocator allocator_;
  std::string name_;

  VarmapT varmap_;

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
        return std::visit(
            [&](auto type) {
              using SelfT = typename TraitT::template VariantType<decltype(type)>;
              if constexpr (!std::is_same_v<void, SelfT>) {
                return sizeof(SelfT) * count;
              }
              return static_cast<std::size_t>(0);
            },
            identifier.config);
      };
    } else {
      using SelfT = typename TraitT::Type;
      m.bytes = sizeof(SelfT) * count;
      m.bytesLayer = [count](const LayerIdentifier& /*identifier*/) {
        return sizeof(SelfT) * count;
      };
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

  void setName(const std::string& name) { this->name_ = name; }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& container : memoryContainer_) {
      container.synchronizeTo(place, stream);
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
    assert(memoryContainer.size() > index);
    return memoryContainer_[index].get(place);
  }

  template <typename HandleT>
  typename HandleT::Type* var(const HandleT& handle,
                              AllocationPlace place = AllocationPlace::Host) {
    return static_cast<typename HandleT::Type*>(varUntyped(varmap_.index(handle), place));
  }

  template <typename StorageT>
  typename StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer.size() > index);
    return static_cast<typename StorageT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename HandleT>
  const typename HandleT::Type* var(const HandleT& handle,
                                    AllocationPlace place = AllocationPlace::Host) const {
    return static_cast<typename HandleT::Type*>(varUntyped(varmap_.index(handle), place));
  }

  template <typename StorageT>
  const typename StorageT::Type* var(AllocationPlace place = AllocationPlace::Host) const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer.size() > index);
    return static_cast<typename StorageT::Type*>(memoryContainer_[index].get(place));
  }

  template <typename StorageT>
  void varSynchronizeTo(AllocationPlace place, void* stream) {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryContainer.size() > index);
    memoryContainer_[index].synchronizeTo(place, stream);
  }

  template <typename StorageT>
  [[nodiscard]] const MemoryInfo& info() const {
    const auto index = varmap_.template index<StorageT>();
    assert(memoryInfo.size() > index);
    return memoryInfo_[index];
  }

  template <typename HandleT>
  [[nodiscard]] const MemoryInfo& info(const HandleT& handle) const {
    const auto index = varmap_.index(handle);
    assert(memoryInfo.size() > index);
    return memoryInfo_[index];
  }

  [[nodiscard]] const MemoryInfo& info(std::size_t index) const {
    assert(memoryInfo.size() > index);
    return memoryInfo_[index];
  }

  auto lookupRef(const StoragePosition& position, AllocationPlace place = AllocationPlace::Host) {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).cellRef(position.cell, place);
  }

  template <typename HandleT>
  auto& lookup(const HandleT& handle,
               const StoragePosition& position,
               AllocationPlace place = AllocationPlace::Host) {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).var(handle, place)[position.cell];
  }

  template <typename StorageT>
  auto& lookup(const StoragePosition& position, AllocationPlace place = AllocationPlace::Host) {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template var<StorageT>(place)[position.cell];
  }

  template <typename HandleT>
  const auto& lookup(const HandleT& handle,
                     const StoragePosition& position,
                     AllocationPlace place = AllocationPlace::Host) const {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).var(handle, place)[position.cell];
  }

  template <typename StorageT>
  const auto& lookup(const StoragePosition& position,
                     AllocationPlace place = AllocationPlace::Host) const {
    assert(position != StoragePosition::NullPosition);
    return layer(position.color).template var<StorageT>(place)[position.cell];
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

  template <typename HandleT>
  void add(HandleT& handle,
           LayerMask mask,
           size_t alignment,
           AllocationMode allocMode,
           bool constant = false,
           std::size_t count = 1) {
    const auto index = varmap_.add(handle);
    addInternal<HandleT>(index, mask, alignment, allocMode, constant, count);
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

    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == MemoryType::Variable) {
        memoryContainer_[var].allocate(allocator_,
                                       memoryInfo_[var].size,
                                       memoryInfo_[var].alignment,
                                       memoryInfo_[var].allocMode);
        memoryContainer_[var].constant = memoryInfo_[var].constant;
      }
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
  }

  // Walks through all leaves, computes the maximum amount of memory for each scratchpad entity,
  // allocates all scratchpads based on evaluated max. scratchpad sizes, and, finally,
  // redistributes scratchpads to all leaves.
  //
  // Note, all scratchpad entities are shared between leaves.
  // Do not update leaves in parallel inside of the same MPI rank while using GPUs.
  void allocateScratchPads() {
    std::vector<std::size_t> sizes(memoryInfo_.size());
    for (auto& leaf : this->leaves()) {
      leaf.findMaxScratchpadSizes(sizes);
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < memoryInfo_.size(); ++var) {
      if (memoryInfo_[var].type == MemoryType::Scratchpad) {
        memoryInfo_[var].size = sizes[var];
        totalSize += sizes[var];
      }
    }
    if (!name_.empty()) {
      logInfo() << "Storage" << name_
                << "; scratchpads:" << UnitByte.formatPrefix(totalSize).c_str();
    }

    for (size_t id = 0; id < memoryInfo_.size(); ++id) {
      if (memoryInfo_[id].type == MemoryType::Scratchpad) {
        memoryContainer_[id].allocate(
            allocator_, memoryInfo_[id].size, memoryInfo_[id].alignment, memoryInfo_[id].allocMode);
      }
    }

    for (auto& leaf : this->leaves()) {
      leaf.setMemoryRegionsForScratchpads(memoryContainer_);
    }
  }

  void touchVariables() {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      for (auto& leaf : this->leaves()) {
        leaf.touchVariables();
      }
    }
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
