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

#include "Lut.hpp"

seissol::initializers::Lut::Lut()
  : m_ltsToMesh(NULL), m_meshToLts(NULL), m_ltsTree(NULL), duplicatedMeshIds(NULL), numberOfDuplicatedMeshIds(0)
{
  for (unsigned i = 0; i < (1 << NUMBER_OF_LAYERS); ++i) {
    m_meshToCell[i] = NULL;
  }
}

seissol::initializers::Lut::~Lut()
{
  delete[] duplicatedMeshIds;
  for (unsigned i = 0; i < (1 << NUMBER_OF_LAYERS); ++i) {
    delete[] m_meshToCell[i];
  }
  delete[] m_meshToLts;
  delete[] m_ltsToMesh;
}

void seissol::initializers::Lut::createLuts(  LTSTree*        ltsTree,
                                              unsigned*       ltsToMesh,
                                              unsigned        numberOfCells,
                                              unsigned        numberOfMeshIds )
{
  assert(numberOfCells == ltsTree->getNumberOfCells());

  m_ltsToMesh = new unsigned[numberOfCells];
  std::copy(ltsToMesh, ltsToMesh + numberOfCells, m_ltsToMesh);
  m_ltsTree = ltsTree;
  
  // meshToCell luts
  for (unsigned var = 0; var < m_ltsTree->getNumberOfVariables(); ++var) {
    LayerMask mask = m_ltsTree->info(var).mask;
    unsigned*& meshToCell = m_meshToCell[mask.to_ulong()];
    if (meshToCell == NULL) {
      unsigned offset = 0;
      meshToCell = new unsigned[numberOfMeshIds];
      std::fill(meshToCell, meshToCell + numberOfMeshIds, std::numeric_limits<unsigned>::max());
      
      for (LTSTree::leaf_iterator it = m_ltsTree->beginLeaf(); it != m_ltsTree->endLeaf(); ++it) {
        if (!it->isMasked(mask)) {
          for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
            unsigned meshId = m_ltsToMesh[it->getLtsIdStart() + cell];
            if (meshId != std::numeric_limits<unsigned>::max()) {
              meshToCell[meshId] = offset + cell;
            }
          }
          offset += it->getNumberOfCells();
        }
      }
    }
  }
  
  // duplicates     
  m_meshToLts               = new unsigned[numberOfMeshIds][MaxDuplicates];
  for (unsigned meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
      m_meshToLts[meshId][dup] = std::numeric_limits<unsigned>::max();
    }    
  }

  unsigned*  numDuplicates  = new unsigned[numberOfMeshIds];
  memset(numDuplicates, 0, numberOfMeshIds * sizeof(unsigned));

  for (unsigned ltsId = 0; ltsId < numberOfCells; ++ltsId) {
    unsigned meshId = m_ltsToMesh[ltsId];
    if (meshId != std::numeric_limits<unsigned>::max()) {
      assert( numDuplicates[meshId] < MaxDuplicates);
      m_meshToLts[meshId][ numDuplicates[meshId]++ ] = ltsId;
    }
  }
  
  numberOfDuplicatedMeshIds = 0;
  for (unsigned meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    if (numDuplicates[meshId] > 1) {
      ++numberOfDuplicatedMeshIds;
    }
  }
  
  duplicatedMeshIds = new unsigned[numberOfDuplicatedMeshIds];
  
  unsigned dupId = 0;
  for (unsigned meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    if (numDuplicates[meshId] > 1) {
      duplicatedMeshIds[dupId++] = meshId;
    }
  }
  
  delete[] numDuplicates;
}
