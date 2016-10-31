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
#include <bitset>
#include <limits>
#include <cstring>

enum LayerType {
  Ghost    = (1 << 0),
  Copy     = (1 << 1),
  Interior = (1 << 2),
  NUMBER_OF_LAYERS
};

namespace seissol {
  namespace initializers {
    typedef std::bitset<NUMBER_OF_LAYERS> LayerMask;

    template<typename T> class Variable;
    class Bucket;
    struct MemoryInfo;
    class Layer;
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
    assert(m_buckets != NULL && m_buckets[handle.index] != NULL);
    return m_buckets[handle.index];
  }
  
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
  
  inline void setBucketSize(Bucket const& handle, size_t size) {
    assert(m_bucketSizes != NULL);
    m_bucketSizes[handle.index] = size;
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
  
  void touchVariables(std::vector<MemoryInfo> const& vars) {
    for (unsigned var = 0; var < vars.size(); ++var) {
      if (!isMasked(vars[var].mask)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (unsigned cell = 0; cell < m_numberOfCells; ++cell) {
          memset(static_cast<char*>(m_vars[var]) + cell * vars[var].bytes, 0, vars[var].bytes);
        }
      }
    }
  }
};

#endif
