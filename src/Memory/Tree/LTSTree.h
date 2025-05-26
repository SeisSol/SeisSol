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
#include "TimeCluster.h"

#include "Memory/MemoryAllocator.h"

#include "Monitoring/Unit.h"
#include "utils/logger.h"

namespace seissol::initializer {

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
    m.bytes = sizeof(typename TraitT::Type);
    m.alignment = alignment;
    m.mask = mask;
    m.allocMode = allocMode;
    m.constant = constant;
    m.type = TraitT::Storage;
    m.index = memoryInfo.size();

    m.filterLayer = [mask](const LayerIdentifier& identifier) {
      return (mask.to_ulong() & identifier.halo) != 0;
    };
    /*
    m.bytesLayer = [](const LayerIdentifier& identifier) {
      return std::visit([&](auto type) {
        return sizeof(StorageT::template TypeVariant<decltype(type)>);
      }, identifier.config);
    };
    */
    m.bytesLayer = [](const LayerIdentifier& identifier) { return sizeof(typename TraitT::Type); };
    memoryInfo.push_back(m);
  }

  public:
  LTSTree() = default;

  ~LTSTree() override = default;

  void setName(const std::string& name) { this->name = name; }

  void synchronizeTo(AllocationPlace place, void* stream) {
    for (auto& container : memoryContainer) {
      container.synchronizeTo(place, stream);
    }
  }

  void setNumberOfTimeClusters(unsigned numberOfTimeCluster) {
    setChildren<TimeCluster>(numberOfTimeCluster);
  }

  void fixate() {
    memoryContainer.resize(memoryInfo.size());
    setPostOrderPointers();
    for (auto& leaf : leaves()) {
      leaf.fixPointers(memoryInfo, typemap, handlemap);
    }
  }

  TimeCluster& child(unsigned index) {
    return *dynamic_cast<TimeCluster*>(m_children[index].get());
  }

  [[nodiscard]] const TimeCluster& child(unsigned index) const {
    return *dynamic_cast<TimeCluster*>(m_children[index].get());
  }

  void* varUntyped(std::size_t index, AllocationPlace place = AllocationPlace::Host) {
    assert(index != std::numeric_limits<unsigned>::max());
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

  [[nodiscard]] const MemoryInfo& info(unsigned index) const {
    assert(memoryInfo.size() > index);
    return memoryInfo[index];
  }

  [[nodiscard]] unsigned getNumberOfVariables() const { return memoryInfo.size(); }

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
    for (unsigned var = 0; var < memoryInfo.size(); ++var) {
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
      memoryContainer[id].allocate(
          allocator, memoryInfo[id].size, memoryInfo[id].alignment, memoryInfo[id].allocMode);
    }

    for (auto& leaf : leaves()) {
      leaf.setMemoryRegionsForScratchpads(memoryContainer);
    }
  }

  void touchVariables() {
    for (auto& leaf : leaves()) {
      leaf.touchVariables(memoryInfo);
    }
  }

  [[nodiscard]] size_t getMaxClusterSize(LayerMask mask = LayerMask()) const {
    size_t maxClusterSize{0};
    for (const auto& leaf : leaves(mask)) {
      const auto currClusterSize = static_cast<size_t>(leaf.getNumberOfCells());
      maxClusterSize = std::max(currClusterSize, maxClusterSize);
    }
    return maxClusterSize;
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_TREE_LTSTREE_H_
