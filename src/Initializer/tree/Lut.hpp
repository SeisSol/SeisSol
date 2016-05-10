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
 * Handles mapping between mesh and cells.
 **/
 
#ifndef INITIALIZER_TREE_LUT_HPP_
#define INITIALIZER_TREE_LUT_HPP_

#include "LTSTree.hpp"

#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace initializers {
    struct DuplicatedCells;
    class Lut;
  }
}

struct seissol::initializers::DuplicatedCells {
  unsigned*  meshIds;
  unsigned (*duplicates)[4]; // contains lts ids
  unsigned   numberOfDuplicates;
  
  unsigned findDuplicateId(unsigned meshId) const {
    // interval [left, right)
    unsigned left = 0;
    unsigned right = numberOfDuplicates;
    while (left < right) {
      unsigned middle = (left+right) / 2;
      if (meshId > meshIds[middle]) {
        left = middle + 1;
      } else if (meshId < meshIds[middle]) {
        right = middle;
      } else {
        return middle;
      }
    }
    return std::numeric_limits<unsigned>::max();
  }
  
  DuplicatedCells() : meshIds(NULL), duplicates(NULL) {}
};

class seissol::initializers::Lut {
private:
  unsigned*                         m_ltsToMesh;
  unsigned*                         m_meshToCell[1 << NUMBER_OF_LAYERS];
  LTSTree*                          m_ltsTree;
  seissol::memory::ManagedAllocator m_allocator;

public:
  DuplicatedCells                   duplicatedCells;
  
  Lut() : m_ltsToMesh(NULL), m_ltsTree(NULL) {
    for (unsigned i = 0; i < (1 << NUMBER_OF_LAYERS); ++i) {
      m_meshToCell[i] = NULL;
    }
  }

  void createLuts(  LTSTree*        ltsTree,
                    unsigned*       ltsToMesh,
                    unsigned        numberOfCells,
                    unsigned        numberOfMeshIds );
  
  inline unsigned meshId(unsigned ltsId) const {
    return m_ltsToMesh[ltsId];
  }
  
  template<typename T>
  unsigned* getMeshToCellLut(Variable<T> const& handle) {
    return m_meshToCell[ m_ltsTree->info(handle.index).mask.to_ulong() ];
  }
  
  template<typename T>
  T& lookup(Variable<T> const& handle, unsigned meshId) {
    unsigned offset = getMeshToCellLut(handle)[meshId];
    return m_ltsTree->var(handle)[offset];
  }
};

#endif
