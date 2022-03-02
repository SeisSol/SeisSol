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

#include "LTSInternalNode.hpp"
#include "TimeCluster.hpp"

#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace initializers {
    class LTSTree;
  }
}

class seissol::initializers::LTSTree : public seissol::initializers::LTSInternalNode {
private:
  void** m_vars;
  void** m_buckets;
  std::vector<MemoryInfo> varInfo;
  std::vector<MemoryInfo> bucketInfo;
  seissol::memory::ManagedAllocator m_allocator;
  std::vector<size_t> variableSizes{};  /*!< sizes of variables within the entire tree in bytes */
  std::vector<size_t> bucketSizes{};    /*!< sizes of buckets within the entire tree in bytes */

#ifdef ACL_DEVICE
  std::vector<MemoryInfo> scratchpadMemInfo{};
  std::vector<size_t> scratchpadMemSizes{};  /*!< sizes of variables within the entire tree in bytes */
  void** scratchpadMemories;
  std::vector<int> scratchpadMemIds{};
#endif  // ACL_DEVICE

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
#ifdef ACL_DEVICE
      it->allocateScratchpadArrays(scratchpadMemInfo.size());
#endif
    }
  }
  
  inline TimeCluster& child(unsigned index) {
    return *static_cast<TimeCluster*>(m_children[index]);
  }
  
  inline TimeCluster const& child(unsigned index) const {
    return *static_cast<TimeCluster*>(m_children[index]);
  }

  template<typename T>
  T* var(Variable<T> const& handle) {
    assert(handle.index != std::numeric_limits<unsigned>::max());
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
    handle.mask = mask;
    MemoryInfo m;
    m.bytes = sizeof(T)*handle.count;
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

#ifdef ACL_DEVICE
  void addScratchpadMemory(ScratchpadMemory& handle, size_t alignment, seissol::memory::Memkind memkind) {
    handle.index = scratchpadMemInfo.size();
    MemoryInfo memoryInfo;
    memoryInfo.alignment = alignment;
    memoryInfo.memkind = memkind;
    scratchpadMemInfo.push_back(memoryInfo);
  }
#endif // ACL_DEVICE
  
  void allocateVariables() {
    m_vars = new void*[varInfo.size()];
    variableSizes.resize(varInfo.size(), 0);

    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->addVariableSizes(varInfo, variableSizes);
    }

    for (unsigned var = 0; var < varInfo.size(); ++var) {
      m_vars[var] = m_allocator.allocateMemory(variableSizes[var], varInfo[var].alignment, varInfo[var].memkind);
    }
    
    std::fill(variableSizes.begin(), variableSizes.end(), 0);
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForVariables(varInfo, m_vars, variableSizes);
      it->addVariableSizes(varInfo, variableSizes);
    }
  }
  
  void allocateBuckets() {
    m_buckets = new void*[bucketInfo.size()];
    bucketSizes.resize(bucketInfo.size(), 0);
    
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->addBucketSizes(bucketSizes);
    }
    
    for (unsigned bucket = 0; bucket < bucketInfo.size(); ++bucket) {
      m_buckets[bucket] = m_allocator.allocateMemory(bucketSizes[bucket], bucketInfo[bucket].alignment, bucketInfo[bucket].memkind);
    }
    
    std::fill(bucketSizes.begin(), bucketSizes.end(), 0);
      for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForBuckets(m_buckets, bucketSizes);
      it->addBucketSizes(bucketSizes);
    }
  }

#ifdef ACL_DEVICE
  // Walks through all leaves, computes the maximum amount of memory for each scratchpad entity,
  // allocates all scratchpads based on evaluated max. scratchpad sizes, and, finally,
  // redistributes scratchpads to all leaves.
  //
  // Note, all scratchpad entities are shared between leaves.
  // Do not update leaves in parallel inside of the same MPI rank while using GPUs.
  void allocateScratchPads() {
    scratchpadMemories = new void*[scratchpadMemInfo.size()];
    scratchpadMemSizes.resize(scratchpadMemInfo.size(), 0);

    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->findMaxScratchpadSizes(scratchpadMemSizes);
    }

    for (size_t id = 0; id < scratchpadMemSizes.size(); ++id) {
      // TODO {ravil}: check whether the assert makes sense
      //assert((scratchpadMemSizes[id] > 0) && "ERROR: scratchpad mem. size is equal to zero");
      scratchpadMemories[id] = m_allocator.allocateMemory(scratchpadMemSizes[id],
                                                          scratchpadMemInfo[id].alignment,
                                                          scratchpadMemInfo[id].memkind);
    }

    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->setMemoryRegionsForScratchpads(scratchpadMemories, scratchpadMemInfo.size());
    }
  }
#endif
  
  void touchVariables() {
    for (LTSTree::leaf_iterator it = beginLeaf(); it != endLeaf(); ++it) {
      it->touchVariables(varInfo);
    }
  }

  const std::vector<size_t>& getVariableSizes() {
    return variableSizes;
  }

  const std::vector<size_t>& getBucketSizes() {
    return bucketSizes;
  }
};

#endif
