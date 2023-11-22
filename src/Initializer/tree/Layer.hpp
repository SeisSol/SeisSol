/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#ifndef INITIALIZER_TREE_LAYER_HPP_
#define INITIALIZER_TREE_LAYER_HPP_

#include "Node.hpp"
#include <Initializer/MemoryAllocator.h>
#include <Initializer/BatchRecorders/DataTypes/ConditionalTable.hpp>
#include "Initializer/DeviceGraph.h"
#include <bitset>
#include <limits>
#include <cstring>
#include <type_traits>
#include <vector>

enum LayerType {
  Ghost    = (1 << 0),
  Copy     = (1 << 1),
  Interior = (1 << 2),
  NUMBER_OF_LAYERS
};

namespace seissol {
  namespace initializers {
    typedef std::bitset<NUMBER_OF_LAYERS> LayerMask;

    template<typename T> struct Variable;
    struct Bucket;
    struct MemoryInfo;
    class Layer;
#ifdef ACL_DEVICE
    struct ScratchpadMemory;
#endif

    enum class AllocationMode {
      HostOnly,
      HostOnlyHBM,
      HostDeviceUnified,
      HostDeviceSplit,
      DeviceOnly
    };

    enum class AllocationPlace {
      Host,
      Device
    };

    struct DualMemoryContainer {
      void* host = nullptr;
      void* device = nullptr;
      AllocationMode allocationMode;
      std::size_t allocationSize;
      std::size_t allocationAlignment;

      void* get(AllocationPlace place) const {
        if (place == AllocationPlace::Host) {
          return host;
        }
        else if (place == AllocationPlace::Device) {
          return device;
        }
        else {
          return nullptr; // should not happen
        }
      }

      void allocate(seissol::memory::ManagedAllocator& allocator, std::size_t size, std::size_t alignment, AllocationMode mode) {
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
          device = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceGlobalMemory);
        }
        if (mode == AllocationMode::HostDeviceUnified) {
          host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceUnifiedMemory);
          device = host;
        }
        if (mode == AllocationMode::HostDeviceSplit) {
          host = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::PinnedMemory);
          device = allocator.allocateMemory(size, alignment, seissol::memory::Memkind::DeviceGlobalMemory);
        }
        allocationMode = mode;
        allocationSize = size;
      }

      void synchronizeTo(AllocationPlace place, void* stream) {
#ifdef ACL_DEVICE
        if (allocationMode == AllocationMode::HostDeviceSplit) {
          if (place == AllocationPlace::Host) {
            device::DeviceInstance::getInstance().api->copyFromAsync(host, device, allocationSize, stream);
          }
          else {
            device::DeviceInstance::getInstance().api->copyToAsync(device, host, allocationSize, stream);
          }

        }
        if (allocationMode == AllocationMode::HostDeviceUnified) {
          if (place == AllocationPlace::Host) {
            device::DeviceInstance::getInstance().api->prefetchUnifiedMemTo(device::Destination::Host, host, allocationSize, stream);
          }
          else {
            device::DeviceInstance::getInstance().api->prefetchUnifiedMemTo(device::Destination::CurrentDevice, device, allocationSize, stream);
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
  }
}

template<typename T>
struct seissol::initializers::Variable {
  unsigned index;
  LayerMask mask;
  unsigned count;
  Variable() : index(std::numeric_limits<unsigned>::max()), count(1) {}
};

struct seissol::initializers::Bucket {
  unsigned index;

  Bucket() : index(std::numeric_limits<unsigned>::max()) {}
};

#ifdef ACL_DEVICE
struct seissol::initializers::ScratchpadMemory : public seissol::initializers::Bucket{};
#endif

struct seissol::initializers::MemoryInfo {
  size_t bytes;
  size_t alignment;
  LayerMask mask;
  // seissol::memory::Memkind memkind;
  AllocationMode allocMode;
};

class seissol::initializers::Layer : public seissol::initializers::Node {
private:
  enum LayerType m_layerType;
  unsigned m_numberOfCells;
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
  Layer() : m_numberOfCells(0) {}
  ~Layer() {  }
  
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

  template<typename T>
  T* var(Variable<T> const& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_vars.size() > handle.index);
    return static_cast<T*>(m_vars[handle.index].get(place));
  }

  void* bucket(Bucket const& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_buckets.size() > handle.index);
    return m_buckets[handle.index].get(place);
  }

#ifdef ACL_DEVICE
  void* getScratchpadMemory(ScratchpadMemory const& handle, AllocationPlace place = AllocationPlace::Host) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_scratchpads.size() > handle.index);
    return (m_scratchpads[handle.index].get(place));
  }
#endif
  
  /// i-th bit of layerMask shall be set if data is masked on the i-th layer
  inline bool isMasked(LayerMask layerMask) const {
    return (LayerMask(m_layerType) & layerMask).any();
  }

  inline void setLayerType(enum LayerType layerType) {
    m_layerType = layerType;
  }
  
  inline enum LayerType getLayerType() const {
    return m_layerType;
  }
  
  inline unsigned getNumberOfCells() const {
    return m_numberOfCells;
  }
  
  inline void setNumberOfCells(unsigned numberOfCells) {
    m_numberOfCells = numberOfCells;
  }
  
  inline void allocatePointerArrays(unsigned numVars, unsigned numBuckets) {
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
  
  inline void setBucketSize(Bucket const& handle, size_t size) {
    assert(m_bucketSizes.size() > handle.index);
    m_bucketSizes[handle.index] = size;
  }

#ifdef ACL_DEVICE
  inline void setScratchpadSize(ScratchpadMemory const& handle, size_t size) {
    assert(m_scratchpadSizes.size() > handle.index);
    m_scratchpadSizes[handle.index] = size;
  }
#endif

  inline size_t getBucketSize(Bucket const& handle) {
    assert(m_bucketSizes.size() > handle.index);
    return m_bucketSizes[handle.index];
    }
  
  void addVariableSizes(std::vector<MemoryInfo> const& vars, std::vector<size_t>& bytes) {
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

  void setMemoryRegionsForVariables(std::vector<MemoryInfo> const& vars, const std::vector<DualMemoryContainer>& memory, const std::vector<size_t>& offsets) {
    assert(m_vars.size() >= vars.size());
    for (unsigned var = 0; var < vars.size(); ++var) {
      if (!isMasked(vars[var].mask)) {
        m_vars[var].offsetFrom(memory[var], offsets[var], m_numberOfCells * vars[var].bytes);
      }
    }
  }

  void setMemoryRegionsForBuckets(const std::vector<DualMemoryContainer>& memory, const std::vector<size_t>& offsets) {
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
  
  void touchVariables(std::vector<MemoryInfo> const& vars) {
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
  template<typename InnerKeyType>
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

  template<typename InnerKeyType>
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
    }
    else {
      return device::DeviceGraphHandle();
    }
  }

  void updateDeviceComputeGraphHandle(GraphKey graphKey,
                                      device::DeviceGraphHandle graphHandle) {
    assert(m_computeGraphHandles.find(graphKey) == m_computeGraphHandles.end()
           && "an entry of hash table must be empty on write");
    if (graphHandle.isInitialized()) {
      m_computeGraphHandles[graphKey] = graphHandle;
    }
  }
#endif // ACL_DEVICE
};

#endif
