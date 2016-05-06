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
 
#ifndef INITIALIZER_TREE_LUT_HPP_
#define INITIALIZER_TREE_LUT_HPP_

#include "LTSTree.hpp"

#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace initializers {
    class Lut;
  }
}

class seissol::initializers::Lut {
private:
  unsigned*                         m_ltsToMesh;
  unsigned*                         m_meshToCell[1 << NUMBER_OF_LAYERS];
  LTSTree*                          m_ltsTree;
  seissol::memory::ManagedAllocator m_allocator;

public:
  Lut() : m_ltsToMesh(NULL), m_ltsTree(NULL) {
    for (unsigned i = 0; i < (1 << NUMBER_OF_LAYERS); ++i) {
      m_meshToCell[i] = NULL;
    }
  }

  void createLuts(  LTSTree*        ltsTree,
                    unsigned*       ltsToMesh,
                    unsigned        numberOfCells,
                    unsigned        numberOfMeshIds ) {
    assert(numberOfCells == ltsTree->getNumberOfCells());

    m_ltsToMesh = static_cast<unsigned*>(m_allocator.allocateMemory(numberOfCells * sizeof(unsigned)));
    std::copy(ltsToMesh, ltsToMesh + numberOfCells, m_ltsToMesh);
    m_ltsTree = ltsTree;
    
    for (unsigned var = 0; var < m_ltsTree->getNumberOfVariables(); ++var) {
      LayerMask mask = m_ltsTree->info(var).mask;
      unsigned*& meshToCell = m_meshToCell[mask.to_ulong()];
      if (meshToCell == NULL) {
        unsigned startLtsId = 0;
        unsigned offset = 0;
        meshToCell = static_cast<unsigned*>(m_allocator.allocateMemory(numberOfMeshIds * sizeof(unsigned)));
        std::fill(meshToCell, meshToCell + numberOfMeshIds, std::numeric_limits<unsigned>::max());
        
        for (LTSTree::leaf_iterator it = m_ltsTree->beginLeaf(); it != m_ltsTree->endLeaf(); ++it) {
          if (!it->isMasked(mask)) {
            for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {      
              unsigned meshId = m_ltsToMesh[startLtsId + cell];
              if (meshId != std::numeric_limits<unsigned>::max()) {
                meshToCell[meshId] = offset + cell;
              }
            }
            offset += it->getNumberOfCells();
          }
          startLtsId += it->getNumberOfCells();
        }
      }
    }
  }
  
  inline unsigned meshId(unsigned ltsId) const {
    return m_ltsToMesh[ltsId];
  }
  
  template<typename T>
  unsigned* getMeshToCellLut(Variable<T>& handle) {
    return m_meshToCell[ m_ltsTree->info(handle.index).mask.to_ulong() ];
  }
  
  template<typename T>
  T& lookup(Variable<T>& handle, unsigned meshId) {
    unsigned offset = getMeshToCellLut(handle)[meshId];
    return m_ltsTree->var(handle)[offset];
  }
};

#endif
