// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#include "Lut.h"
#include <Initializer/Tree/LTSTree.h>
#include <Initializer/Tree/Layer.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

seissol::initializer::Lut::LutsForMask::LutsForMask()
    : ltsToMesh(NULL), duplicatedMeshIds(NULL), numberOfDuplicatedMeshIds(0) {
  for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
    meshToLts[dup] = NULL;
  }
}
seissol::initializer::Lut::LutsForMask::~LutsForMask() {
  delete[] duplicatedMeshIds;
  for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
    delete[] meshToLts[dup];
  }
  delete[] ltsToMesh;
}

void seissol::initializer::Lut::LutsForMask::createLut(LayerMask mask,
                                                       LTSTree* ltsTree,
                                                       const unsigned* globalLtsToMesh,
                                                       unsigned numberOfMeshIds) {
  const unsigned numberOfLtsIds = ltsTree->getNumberOfCells(mask);
  ltsToMesh = new unsigned[numberOfLtsIds];

  // ltsToMesh
  unsigned globalLtsId = 0;
  unsigned offset = 0;
  for (auto it = ltsTree->beginLeaf(); it != ltsTree->endLeaf(); ++it) {
    if (!it->isMasked(mask)) {
      for (unsigned cell = 0; cell < it->getNumberOfCells(); ++cell) {
        const unsigned meshId = globalLtsToMesh[globalLtsId + cell];
        ltsToMesh[offset + cell] = meshId;
      }
      offset += it->getNumberOfCells();
    }
    globalLtsId += it->getNumberOfCells();
  }

  // meshToLts
  for (unsigned dup = 0; dup < MaxDuplicates; ++dup) {
    meshToLts[dup] = new unsigned[numberOfMeshIds];
    std::fill(
        meshToLts[dup], meshToLts[dup] + numberOfMeshIds, std::numeric_limits<unsigned>::max());
  }

  unsigned* numDuplicates = new unsigned[numberOfMeshIds];
  memset(numDuplicates, 0, numberOfMeshIds * sizeof(unsigned));

  for (unsigned ltsId = 0; ltsId < numberOfLtsIds; ++ltsId) {
    const unsigned meshId = ltsToMesh[ltsId];
    if (meshId != std::numeric_limits<unsigned int>::max()) {
      assert(numDuplicates[meshId] < MaxDuplicates);
      meshToLts[numDuplicates[meshId]++][meshId] = ltsId;
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

seissol::initializer::Lut::Lut() : m_ltsTree(NULL), m_meshToClusters(NULL) {}

seissol::initializer::Lut::~Lut() { delete[] m_meshToClusters; }

void seissol::initializer::Lut::createLuts(LTSTree* ltsTree,
                                           unsigned* ltsToMesh,
                                           unsigned numberOfMeshIds) {
  const unsigned numberOfCells = ltsTree->getNumberOfCells();

  m_ltsTree = ltsTree;

  for (unsigned var = 0; var < m_ltsTree->getNumberOfVariables(); ++var) {
    const LayerMask mask = m_ltsTree->info(var).mask;
    LutsForMask& maskedLut = maskedLuts[mask.to_ulong()];
    if (maskedLut.ltsToMesh == NULL) {
      maskedLut.createLut(mask, m_ltsTree, ltsToMesh, numberOfMeshIds);
    }
  }

  struct LayerOffset {
    unsigned offsetGhost;
    unsigned offsetCopy;
    unsigned offsetInterior;
  };

  auto* clusters = new unsigned[m_ltsTree->numChildren() + 1];
  // Store number of cells in layers for each timecluster.
  // Note that the index 0 is the first cluster, unlike as in clusters.
  auto clustersLayerOffset = std::vector<LayerOffset>(m_ltsTree->numChildren());
  clusters[0] = 0;
  for (unsigned tc = 0; tc < m_ltsTree->numChildren(); ++tc) {
    auto& cluster = m_ltsTree->child(tc);
    clusters[tc + 1] = clusters[tc] + cluster.getNumberOfCells();
    // For each cluster, we first store the Ghost cells, then the Copy cells and finally the
    // Interior cells.
    auto offsetGhost = 0U;
    auto offsetCopy = cluster.child<Ghost>().getNumberOfCells();
    auto offsetInterior = offsetCopy + cluster.child<Copy>().getNumberOfCells();
    clustersLayerOffset[tc] = LayerOffset{offsetGhost, offsetCopy, offsetInterior};
  }

  m_meshToClusters = new unsigned[numberOfMeshIds];
  m_meshToLayer.resize(numberOfMeshIds);
  unsigned cluster = 0;
  unsigned curClusterElements = 0;
  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    if (cell >= clusters[cluster + 1]) {
      curClusterElements = 0;
      ++cluster;
      assert(cluster <= ltsTree->numChildren());
    } else {
      ++curClusterElements;
    }
    if (ltsToMesh[cell] != std::numeric_limits<unsigned>::max()) {
      const auto meshId = ltsToMesh[cell];
      m_meshToClusters[meshId] = cluster;
      const auto& layerOffsets = clustersLayerOffset[cluster];
      if (curClusterElements >= layerOffsets.offsetInterior) {
        m_meshToLayer[meshId] = Interior;
      } else if (curClusterElements >= layerOffsets.offsetCopy) {
        m_meshToLayer[meshId] = Copy;
      } else if (curClusterElements >= layerOffsets.offsetGhost) {
        m_meshToLayer[meshId] = Ghost;
      } else {
        throw std::logic_error("Can't tell which layer the meshid belongs.");
      }
    }
  }

  delete[] clusters;
}
