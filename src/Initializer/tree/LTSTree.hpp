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
 * Tree for managing lts data.
 **/
 
#ifndef INITIALIZER_TREE_LTSTREE_HPP_
#define INITIALIZER_TREE_LTSTREE_HPP_

#include "Node.hpp"
#include "Layer.hpp"
#include "TimeCluster.hpp"

#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace initializers {
    class LTSTree;
  }
}

class seissol::initializers::LTSTree : public seissol::initializers::Node {
private:
  void** m_vars;
  void** m_buckets;
  std::vector<MemoryInfo> varInfo;
  std::vector<MemoryInfo> bucketInfo;
  seissol::memory::ManagedAllocator m_allocator;

public:
  LTSTree() : m_vars(NULL), m_buckets(NULL) {}
  
  ~LTSTree() { delete[] m_vars; delete[] m_buckets; }
  
  void setNumberOfTimeClusters(unsigned numberOfTimeCluster) {
    setChildren<TimeCluster>(numberOfTimeCluster);
  }
  
  void fixate() {
    setPostOrderPointers();
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->allocatePointerArrays(varInfo.size(), bucketInfo.size());
    }
  }
  
  inline TimeCluster& child(unsigned index) {
    return *static_cast<TimeCluster*>(m_children[index]);
  }
  
  /// \todo remove?
  template<typename T>
  T* var(Variable<T> const& handle) {
    assert(m_vars != NULL/* && m_vars[handle.index] != NULL*/);
    return static_cast<T*>(m_vars[handle.index]);
  }

  MemoryInfo const& info(unsigned index) const {
    return varInfo[index];
  }
  
  inline unsigned getNumberOfVariables() const {
    return varInfo.size();
  }
  
  template<typename T>
  void addVar(Variable<T>& handle, LayerMask mask, size_t alignment, seissol::memory::Memkind memkind) {
    handle.index = varInfo.size();
    MemoryInfo m;
    m.bytes = sizeof(T);
    m.alignment = alignment;
    m.mask = mask;
    m.memkind = memkind;
    varInfo.push_back(m);
  }
  
  void addBucket(Bucket& handle, size_t alignment, seissol::memory::Memkind memkind) {
    handle.index = bucketInfo.size();
    MemoryInfo m;
    m.alignment = alignment;
    m.memkind = memkind;
    bucketInfo.push_back(m);
  }
  
  void allocateMemory() {
    m_vars = new void*[varInfo.size()];
    m_buckets = new void*[bucketInfo.size()];
    std::vector<size_t> variableSizes(varInfo.size(), 0);
    std::vector<size_t> bucketSizes(bucketInfo.size(), 0);
    
    unsigned ltsIdStart = 0;
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->addVariableSizes(varInfo, variableSizes);
      it->addBucketSizes(bucketSizes);
      it->setLtsIdStart(ltsIdStart);
      ltsIdStart += it->getNumberOfCells();
    }

    for (unsigned var = 0; var < varInfo.size(); ++var) {
      m_vars[var] = m_allocator.allocateMemory(variableSizes[var], varInfo[var].alignment, varInfo[var].memkind);
    }
    for (unsigned bucket = 0; bucket < bucketInfo.size(); ++bucket) {
      m_buckets[bucket] = m_allocator.allocateMemory(bucketSizes[bucket], bucketInfo[bucket].alignment, bucketInfo[bucket].memkind);
    }
    
    std::fill(variableSizes.begin(), variableSizes.end(), 0);
    std::fill(bucketSizes.begin(), bucketSizes.end(), 0);
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForVariables(varInfo, m_vars, variableSizes);
      it->addVariableSizes(varInfo, variableSizes);
      it->setMemoryRegionsForBuckets(m_buckets, bucketSizes);
      it->addBucketSizes(bucketSizes);
    }
  }
  
  void touch() {
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->touch(varInfo);
    }
  }
  
  unsigned getNumberOfCells() {
    unsigned numCells = 0;
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      numCells += it->getNumberOfCells();
    }
    return numCells;
  }
  
  class leaf_iterator : public iterator {
    friend class LTSTree;

  private:
    LayerMask m_layerMask;
    
    inline void nextLeaf() {
      do {
        m_node = m_node->m_next;
      } while (m_node != NULL && !m_node->isLeaf());
    }

    // m_node must point to a leaf or NULL
    inline void skipMaskedLayer() {
      while (m_node != NULL && operator*().isMasked(m_layerMask)) {
        nextLeaf();
      }
    }

  public:
    leaf_iterator() : iterator() {}
    leaf_iterator(iterator const& it, LayerMask layerMask) : iterator(it), m_layerMask(layerMask) {}

    inline iterator& operator++() {
      nextLeaf();
      skipMaskedLayer();
      return *this;
    }
    
    inline Layer& operator*() {
      return *static_cast<Layer*>(m_node);
    }
    
    inline Layer* operator->() {
      return static_cast<Layer*>(m_node);
    }
  };

  inline leaf_iterator beginLeaf(LayerMask layerMask = LayerMask()) {
    leaf_iterator it = leaf_iterator(begin(), layerMask);
    it.skipMaskedLayer();
    return it;
  }
  
  inline leaf_iterator endLeaf() {
    return leaf_iterator();
  }
};

#endif
