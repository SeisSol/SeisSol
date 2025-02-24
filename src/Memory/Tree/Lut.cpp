// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Lut.h"
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <cassert>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

void seissol::initializer::Lut::LutsForMask::createLut(LayerMask mask,
                                                       LTSTree* ltsTree,
                                                       const unsigned* globalLtsToMesh,
                                                       unsigned numberOfMeshIds) {
  const unsigned numberOfLtsIds = ltsTree->getNumberOfCells(mask);
  ltsToMesh.resize(numberOfLtsIds);

  // ltsToMesh
  unsigned globalLtsId = 0;
  unsigned offset = 0;
  for (const auto& layer : ltsTree->leaves()) {
    if (!layer.isMasked(mask)) {
      for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
        const unsigned meshId = globalLtsToMesh[globalLtsId + cell];
        ltsToMesh[offset + cell] = meshId;
      }
      offset += layer.getNumberOfCells();
    }
    globalLtsId += layer.getNumberOfCells();
  }

  // meshToLts
  for (auto& meshToLt : meshToLts) {
    meshToLt.resize(numberOfMeshIds, std::numeric_limits<unsigned>::max());
  }

  std::vector<unsigned> numDuplicates(numberOfMeshIds);

  for (unsigned ltsId = 0; ltsId < numberOfLtsIds; ++ltsId) {
    const unsigned meshId = ltsToMesh[ltsId];
    if (meshId != std::numeric_limits<unsigned int>::max()) {
      assert(numDuplicates[meshId] < MaxDuplicates);
      meshToLts[numDuplicates[meshId]++][meshId] = ltsId;
    }
  }

  std::size_t numberOfDuplicatedMeshIds = 0;
  for (unsigned meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    if (numDuplicates[meshId] > 1) {
      ++numberOfDuplicatedMeshIds;
    }
  }

  duplicatedMeshIds.resize(numberOfDuplicatedMeshIds);

  unsigned dupId = 0;
  for (unsigned meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    if (numDuplicates[meshId] > 1) {
      duplicatedMeshIds[dupId++] = meshId;
    }
  }
}

seissol::initializer::Lut::Lut() = default;

void seissol::initializer::Lut::createLuts(LTSTree* ltsTree,
                                           unsigned* ltsToMesh,
                                           unsigned numberOfMeshIds) {
  const unsigned numberOfCells = ltsTree->getNumberOfCells();

  m_ltsTree = ltsTree;

  for (unsigned var = 0; var < m_ltsTree->getNumberOfVariables(); ++var) {
    const LayerMask mask = m_ltsTree->info(var).mask;
    LutsForMask& maskedLut = maskedLuts[mask.to_ulong()];
    if (maskedLut.ltsToMesh.empty()) {
      maskedLut.createLut(mask, m_ltsTree, ltsToMesh, numberOfMeshIds);
    }
  }

  struct LayerOffset {
    unsigned offsetGhost;
    unsigned offsetCopy;
    unsigned offsetInterior;
  };

  std::vector<unsigned> clusters(m_ltsTree->numChildren() + 1);
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

  m_meshToClusters.resize(numberOfMeshIds);
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
}
