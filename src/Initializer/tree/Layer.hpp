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
#include <bitset>
#include <limits>
#include <cstring>
#include <type_traits>

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
  seissol::memory::Memkind memkind;
};

class seissol::initializers::Layer : public seissol::initializers::Node {
private:
  enum LayerType m_layerType;
  unsigned m_numberOfCells;
  void** m_vars;
  void** m_buckets;
  size_t* m_bucketSizes;

#ifdef ACL_DEVICE
  void** m_scratchpads{};
  size_t* m_scratchpadSizes{};
  ConditionalPointersToRealsTable m_conditionalPointersToRealsTable{};
  DrConditionalPointersToRealsTable m_drConditionalPointersToRealsTable{};
  ConditionalIndicesTable m_conditionalIndicesTable;
#endif

public:
  Layer() : m_numberOfCells(0), m_vars(NULL), m_buckets(NULL), m_bucketSizes(NULL) {}
  ~Layer() { delete[] m_vars; delete[] m_buckets; delete[] m_bucketSizes; }
  
  template<typename T>
  T* var(Variable<T> const& handle) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_vars != NULL/* && m_vars[handle.index] != NULL*/);
    return static_cast<T*>(m_vars[handle.index]);
  }

  void* bucket(Bucket const& handle) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_buckets != nullptr && m_buckets[handle.index] != nullptr);
    return m_buckets[handle.index];
  }

#ifdef ACL_DEVICE
  void* getScratchpadMemory(ScratchpadMemory const& handle) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
    assert(m_scratchpads != NULL/* && m_vars[handle.index] != NULL*/);
    return (m_scratchpads[handle.index]);
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
    assert(m_vars == NULL && m_buckets == NULL && m_bucketSizes == NULL);    
    m_vars = new void*[numVars];
    std::fill(m_vars, m_vars + numVars, static_cast<void*>(NULL));
    m_buckets = new void*[numBuckets];
    std::fill(m_buckets, m_buckets + numBuckets, static_cast<void*>(NULL));
    m_bucketSizes = new size_t[numBuckets];
    std::fill(m_bucketSizes, m_bucketSizes + numBuckets, 0);
  }

#ifdef ACL_DEVICE
  inline void allocateScratchpadArrays(unsigned numScratchPads) {
    assert(m_scratchpads == nullptr && m_scratchpadSizes == nullptr);

    m_scratchpads = new void*[numScratchPads];
    std::fill(m_scratchpads, m_scratchpads + numScratchPads, nullptr);

    m_scratchpadSizes = new size_t[numScratchPads];
    std::fill(m_scratchpadSizes, m_scratchpadSizes + numScratchPads, 0);
  }
#endif
  
  inline void setBucketSize(Bucket const& handle, size_t size) {
    assert(m_bucketSizes != NULL);
    m_bucketSizes[handle.index] = size;
  }

#ifdef ACL_DEVICE
  inline void setScratchpadSize(ScratchpadMemory const& handle, size_t size) {
    assert(m_scratchpadSizes != NULL);
    m_scratchpadSizes[handle.index] = size;
  }
#endif

  inline size_t getBucketSize(Bucket const& handle) {
    assert(m_bucketSizes != nullptr);
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

  void setMemoryRegionsForVariables(std::vector<MemoryInfo> const& vars, void** memory, std::vector<size_t>& offsets) {
    assert(m_vars != NULL);
    for (unsigned var = 0; var < vars.size(); ++var) {
      if (!isMasked(vars[var].mask)) {
        m_vars[var] = static_cast<char*>(memory[var]) + offsets[var];
      }
    }
  }

  void setMemoryRegionsForBuckets(void** memory, std::vector<size_t>& offsets) {
    assert(m_buckets != NULL);
    for (unsigned bucket = 0; bucket < offsets.size(); ++bucket) {
      m_buckets[bucket] = static_cast<char*>(memory[bucket]) + offsets[bucket];
    }
  }

#ifdef ACL_DEVICE
  void setMemoryRegionsForScratchpads(void** memory, size_t numScratchPads) {
    assert(m_scratchpads != NULL);
    for (size_t id = 0; id < numScratchPads; ++id) {
      m_scratchpads[id] = static_cast<char*>(memory[id]);
    }
  }
#endif
  
  void touchVariables(std::vector<MemoryInfo> const& vars) {
    for (unsigned var = 0; var < vars.size(); ++var) {

      // NOTE: we don't touch device global memory because it is in a different address space
      // we will do deep-copy from the host to a device later on
      if (!isMasked(vars[var].mask) && (vars[var].memkind != seissol::memory::DeviceGlobalMemory)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (unsigned cell = 0; cell < m_numberOfCells; ++cell) {
          memset(static_cast<char*>(m_vars[var]) + cell * vars[var].bytes, 0, vars[var].bytes);
        }
      }
    }
  }

#ifdef ACL_DEVICE
  template<typename KeyType>
  auto& getConditionalTable() {
    if constexpr (std::is_same_v<KeyType, inner_keys::Wp>) {
      return m_conditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<KeyType, inner_keys::Dr>) {
      return m_drConditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<KeyType, inner_keys::Indices>) {
      return m_conditionalIndicesTable;
    }
  }

  template<typename KeyType>
  const auto& getConditionalTable() const {
    if constexpr (std::is_same_v<KeyType, inner_keys::Wp>) {
      return m_conditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<KeyType, inner_keys::Dr>) {
      return m_drConditionalPointersToRealsTable;
    }

    if constexpr (std::is_same_v<KeyType, inner_keys::Indices>) {
      return m_conditionalIndicesTable;
    }
  }
#endif // ACL_DEVICE
};

#endif
