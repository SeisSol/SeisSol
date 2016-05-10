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

namespace seissol {
  namespace initializers {
    class Lut;
  }
}

class seissol::initializers::Lut {
public:
  static unsigned const             MaxDuplicates = 4;

private:
  unsigned*                         m_ltsToMesh;
  unsigned                        (*m_meshToLts)[MaxDuplicates];
  unsigned*                         m_meshToCell[1 << NUMBER_OF_LAYERS];
  LTSTree*                          m_ltsTree;

public:
  unsigned*                         duplicatedMeshIds;
  unsigned                          numberOfDuplicatedMeshIds;
  
  Lut();  
  ~Lut();

  void createLuts(  LTSTree*        ltsTree,
                    unsigned*       ltsToMesh,
                    unsigned        numberOfCells,
                    unsigned        numberOfMeshIds );
  
  inline unsigned meshId(unsigned ltsId) const {
    return m_ltsToMesh[ltsId];
  }
  
  inline unsigned (&ltsIds(unsigned meshId) const)[MaxDuplicates] {
    return m_meshToLts[meshId];
  }
  
  template<typename T>
  unsigned const* getMeshToCellLut(Variable<T> const& handle) const {
    return m_meshToCell[ m_ltsTree->info(handle.index).mask.to_ulong() ];
  }
  
  template<typename T>
  T& lookup(Variable<T> const& handle, unsigned meshId) const {
    unsigned offset = getMeshToCellLut(handle)[meshId];
    return m_ltsTree->var(handle)[offset];
  }
};

#endif
