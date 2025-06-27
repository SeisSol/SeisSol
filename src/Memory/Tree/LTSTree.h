// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_
#define SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_

#include "LTSInternalNode.h"
#include "Layer.h"

#include "Memory/MemoryAllocator.h"

#include "Monitoring/Unit.h"
#include "utils/logger.h"
#include <Memory/Tree/Colormap.h>
#include <type_traits>

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

struct LayerDefinition {
  LayerIdentifier identifier;
  std::size_t size;
};

class LTSTree : public LTSInternalNode {
  private:
  std::vector<DualMemoryContainer> memoryContainer;
  std::vector<MemoryInfo> memoryInfo;

  seissol::memory::ManagedAllocator allocator;
  std::string name;

  std::unordered_map<std::type_index, std::size_t> typemap;
  std::unordered_map<int*, std::size_t> handlemap;

  template <typename TraitT>
  void addInternal(LayerMask mask,
                   size_t alignment,
                   AllocationMode allocMode,
                   bool constant,
                   std::size_t count) {
    MemoryInfo m;
    m.alignment = alignment;
    m.mask = mask;
    m.allocMode = allocMode;
    m.constant = constant;
    m.type = TraitT::Storage;
    m.index = memoryInfo.size();

    if constexpr (TraitT::IsVariant) {
      m.bytes = 0;
      m.bytesLayer = [](const LayerIdentifier& identifier) {
        return std::visit(
            [&](auto type) {
              using SelfT = typename TraitT::template VariantType<decltype(type)>;
              if constexpr (!std::is_same_v<void, SelfT>) {
                return sizeof(SelfT);
              }
              return static_cast<std::size_t>(0);
            },
            identifier.config);
      };
    } else {
      using SelfT = typename TraitT::Type;
      m.bytes = sizeof(SelfT);
      m.bytesLayer = [](const LayerIdentifier& identifier) {
        if constexpr (!std::is_same_v<void, SelfT>) {
          return sizeof(SelfT);
        }
        return static_cast<std::size_t>(0);
      };
    }

    const auto bytesLayer = m.bytesLayer;

    m.filterLayer = [mask, bytesLayer](const LayerIdentifier& identifier) {
      return mask.test(static_cast<int>(identifier.halo)) || bytesLayer(identifier) == 0;
    };

    memoryInfo.push_back(m);
  }

  std::size_t timeClusters;
  std::vector<ConfigVariant> configs;

  std::optional<LTSColorMap> map;

  public:
  LTSTree() = default;

  ~LTSTree() override = default;

  const LTSColorMap& getColorMap() const { return map.value(); }

  void setName(const std::string& name) { this->name = name; }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& container : memoryContainer) {
      container.synchronizeTo(place, stream);
    }
  }

  void initialize(const LTSColorMap& map, const std::vector<std::size_t>& sizes) {
    this->map = map;
    memoryContainer.resize(memoryInfo.size());
    setChildren<Layer>(map.size());
    setPostOrderPointers();
    assert(map.size() == sizes.size());
    for (std::size_t i = 0; i < sizes.size(); ++i) {
      layer(i).fixPointers(sizes[i], i, map.argument(i), memoryInfo, typemap, handlemap);
    }
  }

  std::size_t numTimeClusters() const { return timeClusters; }

  std::vector<ConfigVariant> getConfigs() const { return configs; }

  Layer& layer(std::size_t index) { return *dynamic_cast<Layer*>(m_children[index].get()); }

  [[nodiscard]] const Layer& layer(std::size_t index) const {
    return *dynamic_cast<Layer*>(m_children[index].get());
  }

  Layer& layer(const LayerIdentifier& id) { return layer(map.value().colorId(id)); }

  [[nodiscard]] const Layer& layer(const LayerIdentifier& id) const {
    return layer(map.value().colorId(id));
  }

  void* varUntyped(std::size_t index, AllocationPlace place = AllocationPlace::Host) {
    assert(index != std::numeric_limits<std::size_t>::max());
    assert(memoryContainer.size() > index);
    return memoryContainer[index].get(place);
  }

  template <typename HandleT>
  typename HandleT::Type* var(const HandleT& handle,
                              AllocationPlace place = AllocationPlace::Host) {
    return static_cast<typename HandleT::Type*>(varUntyped(handlemap.at(handle.pointer()), place));
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

  template <typename StorageT>
  [[nodiscard]] const MemoryInfo& info() const {
    const auto index = typemap.at(std::type_index(typeid(StorageT)));
    assert(memoryInfo.size() > index);
    return memoryInfo[index];
  }

  template <typename HandleT>
  [[nodiscard]] const MemoryInfo& info(const HandleT& handle) const {
    const auto index = handlemap.at(handle.pointer());
    assert(memoryInfo.size() > index);
    return memoryInfo[index];
  }

  [[nodiscard]] const MemoryInfo& info(std::size_t index) const {
    assert(memoryInfo.size() > index);
    return memoryInfo[index];
  }

  [[nodiscard]] std::size_t getNumberOfVariables() const { return memoryInfo.size(); }

  template <typename StorageT>
  void add(LayerMask mask,
           size_t alignment,
           AllocationMode allocMode,
           bool constant = false,
           std::size_t count = 1) {
    typemap[std::type_index(typeid(StorageT))] = memoryInfo.size();
    addInternal<StorageT>(mask, alignment, allocMode, constant, count);
  }

  template <typename HandleT>
  void add(HandleT& handle,
           LayerMask mask,
           size_t alignment,
           AllocationMode allocMode,
           bool constant = false,
           std::size_t count = 1) {
    handlemap[handle.pointer()] = memoryInfo.size();
    addInternal<HandleT>(mask, alignment, allocMode, constant, count);
  }

  void allocateTouchVariables() {
    allocateVariables();
    touchVariables();
  }

  void allocateVariables() {
    std::vector<std::size_t> sizes(memoryInfo.size());
    for (auto& leaf : leaves()) {
      leaf.addVariableSizes(memoryInfo, sizes);
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < memoryInfo.size(); ++var) {
      if (memoryInfo[var].type == MemoryType::Variable) {
        memoryInfo[var].size = sizes[var];
        totalSize += sizes[var];
      }
    }
    if (!name.empty()) {
      logInfo() << "Storage" << name << "; variables:" << UnitByte.formatPrefix(totalSize).c_str();
    }

    for (std::size_t var = 0; var < memoryInfo.size(); ++var) {
      if (memoryInfo[var].type == MemoryType::Variable) {
        memoryContainer[var].allocate(
            allocator, memoryInfo[var].size, memoryInfo[var].alignment, memoryInfo[var].allocMode);
        memoryContainer[var].constant = memoryInfo[var].constant;
      }
    }

    std::fill(sizes.begin(), sizes.end(), 0);
    for (auto& leaf : leaves()) {
      leaf.setMemoryRegionsForVariables(memoryInfo, memoryContainer, sizes);
      leaf.addVariableSizes(memoryInfo, sizes);
    }
  }

  void allocateBuckets() {
    std::vector<std::size_t> sizes(memoryInfo.size());
    for (auto& leaf : leaves()) {
      leaf.addBucketSizes(memoryInfo, sizes);
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < memoryInfo.size(); ++var) {
      if (memoryInfo[var].type == MemoryType::Bucket) {
        memoryInfo[var].size = sizes[var];
        totalSize += sizes[var];
      }
    }
    if (!name.empty()) {
      logInfo() << "Storage" << name << "; buckets:" << UnitByte.formatPrefix(totalSize).c_str();
    }

    for (std::size_t bucket = 0; bucket < memoryInfo.size(); ++bucket) {
      if (memoryInfo[bucket].type == MemoryType::Bucket) {
        memoryContainer[bucket].allocate(allocator,
                                         memoryInfo[bucket].size,
                                         memoryInfo[bucket].alignment,
                                         memoryInfo[bucket].allocMode);
      }
    }

    std::fill(sizes.begin(), sizes.end(), 0);
    for (auto& leaf : leaves()) {
      leaf.setMemoryRegionsForBuckets(memoryContainer, sizes);
      leaf.addBucketSizes(memoryInfo, sizes);
    }
  }

  // Walks through all leaves, computes the maximum amount of memory for each scratchpad entity,
  // allocates all scratchpads based on evaluated max. scratchpad sizes, and, finally,
  // redistributes scratchpads to all leaves.
  //
  // Note, all scratchpad entities are shared between leaves.
  // Do not update leaves in parallel inside of the same MPI rank while using GPUs.
  void allocateScratchPads() {
    std::vector<std::size_t> sizes(memoryInfo.size());
    for (auto& leaf : leaves()) {
      leaf.findMaxScratchpadSizes(memoryInfo, sizes);
    }

    std::size_t totalSize = 0;
    for (std::size_t var = 0; var < memoryInfo.size(); ++var) {
      if (memoryInfo[var].type == MemoryType::Scratchpad) {
        memoryInfo[var].size = sizes[var];
        totalSize += sizes[var];
      }
    }
    if (!name.empty()) {
      logInfo() << "Storage" << name
                << "; scratchpads:" << UnitByte.formatPrefix(totalSize).c_str();
    }

    for (size_t id = 0; id < memoryInfo.size(); ++id) {
      if (memoryInfo[id].type == MemoryType::Scratchpad) {
        memoryContainer[id].allocate(
            allocator, memoryInfo[id].size, memoryInfo[id].alignment, memoryInfo[id].allocMode);
      }
    }

    for (auto& leaf : leaves()) {
      leaf.setMemoryRegionsForScratchpads(memoryContainer);
    }
  }

  void touchVariables() {
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      for (auto& leaf : leaves()) {
        leaf.touchVariables(memoryInfo);
      }
    }
  }

  [[nodiscard]] size_t getMaxClusterSize(LayerMask mask = LayerMask()) const {
    size_t maxClusterSize{0};
    for (const auto& leaf : leaves(mask)) {
      const auto currClusterSize = static_cast<size_t>(leaf.size());
      maxClusterSize = std::max(currClusterSize, maxClusterSize);
    }
    return maxClusterSize;
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_
