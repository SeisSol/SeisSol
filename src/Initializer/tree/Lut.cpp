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

seissol::initializers::Lut::LutsForMask::LutsForMask()
  : ltsToMesh(NULL), duplicatedMeshIds(NULL), numberOfDuplicatedMeshIds(0)
{
  for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
    meshToLts[dup] = NULL;
  }
}
seissol::initializers::Lut::LutsForMask::~LutsForMask()
{
  delete[] duplicatedMeshIds;
  for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
    delete[] meshToLts[dup];
  }
  delete[] ltsToMesh;
}

void seissol::initializers::Lut::LutsForMask::createLut(  LayerMask mask,
                                                          LTSTree*  ltsTree,
                                                          unsigned* globalLtsToMesh,
                                                          unsigned  numberOfMeshIds )
{
  unsigned numberOfLtsIds = ltsTree->getNumberOfCells(mask);
  ltsToMesh = new unsigned[numberOfLtsIds];
  
  // ltsToMesh
  unsigned globalLtsId = 0;
  unsigned offset = 0;
  for (LTSTree::leaf_iterator it = ltsTree->beginLeaf(); it != ltsTree->endLeaf(); ++it) {
    if (!it->isMasked(mask)) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        unsigned meshId = globalLtsToMesh[globalLtsId + cell];
        ltsToMesh[offset + cell] = meshId;
      }
      offset += it->getNumberOfCells();
    }
    globalLtsId += it->getNumberOfCells();
  }

  // meshToLts
  for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
    meshToLts[dup] = new unsigned[numberOfMeshIds];
    std::fill(meshToLts[dup], meshToLts[dup] + numberOfMeshIds, std::numeric_limits<unsigned>::max());
  }

  unsigned* numDuplicates = new unsigned[numberOfMeshIds];
  memset(numDuplicates, 0, numberOfMeshIds * sizeof(unsigned));

  for (unsigned ltsId = 0; ltsId < numberOfLtsIds; ++ltsId) {
    unsigned meshId = ltsToMesh[ltsId];
    if (meshId != std::numeric_limits<unsigned int>::max()) {
      assert( numDuplicates[meshId] < MaxDuplicates);
      meshToLts[ numDuplicates[meshId]++ ][meshId] = ltsId;
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

seissol::initializers::Lut::Lut()
  : m_ltsTree(NULL), m_meshToClusters(NULL)
{
}

seissol::initializers::Lut::~Lut()
{
  delete[] m_meshToClusters;
}

void seissol::initializers::Lut::createLuts(  LTSTree*        ltsTree,
                                              unsigned*       ltsToMesh,
                                              unsigned        numberOfMeshIds )
{
  unsigned numberOfCells = ltsTree->getNumberOfCells();

  m_ltsTree = ltsTree;

  for (unsigned var = 0; var < m_ltsTree->getNumberOfVariables(); ++var) {
    LayerMask mask = m_ltsTree->info(var).mask;    
    LutsForMask& maskedLut = maskedLuts[mask.to_ulong()];
    if (maskedLut.ltsToMesh == NULL) {
      maskedLut.createLut(mask, m_ltsTree, ltsToMesh, numberOfMeshIds );
    }
  }
  
  unsigned* clusters = new unsigned[m_ltsTree->numChildren()+1];
  clusters[0] = 0;
  for (unsigned tc = 0; tc < m_ltsTree->numChildren(); ++tc) {
    clusters[tc+1] = clusters[tc] + m_ltsTree->child(tc).getNumberOfCells();
  }
  
  m_meshToClusters = new unsigned[numberOfMeshIds];
  unsigned cluster = 0;
  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    if (cell >= clusters[cluster+1]) {
      ++cluster;
      assert(cluster <= ltsTree->numChildren());
    }
    if (ltsToMesh[cell] != std::numeric_limits<unsigned>::max()) {
      m_meshToClusters[ ltsToMesh[cell] ] = cluster;
    }
  }
  
  delete[] clusters;
}
