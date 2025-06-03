// SPDX-FileCopyrightText: 2016 SeisSol Group
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
                                                       const std::size_t* globalLtsToMesh,
                                                       std::size_t numberOfMeshIds) {
  const std::size_t numberOfLtsIds = ltsTree->size(mask);
  ltsToMesh.resize(numberOfLtsIds);

  // ltsToMesh
  std::size_t globalLtsId = 0;
  std::size_t offset = 0;
  for (const auto& layer : ltsTree->leaves()) {
    // if (!layer.isMasked(mask)) {
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const std::size_t meshId = globalLtsToMesh[globalLtsId + cell];
      ltsToMesh[offset + cell] = meshId;
    }
    offset += layer.size();
    //}
    globalLtsId += layer.size();
  }

  // meshToLts
  for (auto& meshToLt : meshToLts) {
    meshToLt.resize(numberOfMeshIds, std::numeric_limits<std::size_t>::max());
  }

  std::vector<std::size_t> numDuplicates(numberOfMeshIds);

  for (std::size_t ltsId = 0; ltsId < numberOfLtsIds; ++ltsId) {
    const std::size_t meshId = ltsToMesh[ltsId];
    if (meshId != std::numeric_limits<std::size_t>::max()) {
      assert(numDuplicates[meshId] < MaxDuplicates);
      meshToLts[numDuplicates[meshId]++][meshId] = ltsId;
    }
  }

  std::size_t numberOfDuplicatedMeshIds = 0;
  for (std::size_t meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    if (numDuplicates[meshId] > 1) {
      ++numberOfDuplicatedMeshIds;
    }
  }

  duplicatedMeshIds.resize(numberOfDuplicatedMeshIds);

  std::size_t dupId = 0;
  for (std::size_t meshId = 0; meshId < numberOfMeshIds; ++meshId) {
    if (numDuplicates[meshId] > 1) {
      duplicatedMeshIds[dupId++] = meshId;
    }
  }
}

seissol::initializer::Lut::Lut() = default;

void seissol::initializer::Lut::createLuts(LTSTree* ltsTree,
                                           std::size_t* ltsToMesh,
                                           std::size_t numberOfMeshIds) {
  const std::size_t numberOfCells = ltsTree->size();

  m_ltsTree = ltsTree;

  for (std::size_t var = 0; var < m_ltsTree->getNumberOfVariables(); ++var) {
    const LayerMask mask = m_ltsTree->info(var).mask;
    LutsForMask& maskedLut = maskedLuts[mask.to_ulong()];
    if (maskedLut.ltsToMesh.empty()) {
      maskedLut.createLut(mask, m_ltsTree, ltsToMesh, numberOfMeshIds);
    }
  }

  struct LayerOffset {
    std::size_t offsetGhost;
    std::size_t offsetCopy;
    std::size_t offsetInterior;
  };

  std::vector<std::size_t> clusters(m_ltsTree->numTimeClusters() + 1);
  // Store number of cells in layers for each timecluster.
  // Note that the index 0 is the first cluster, unlike as in clusters.
  auto clustersLayerOffset = std::vector<LayerOffset>(m_ltsTree->numTimeClusters());
  clusters[0] = 0;
  for (std::size_t tc = 0; tc < m_ltsTree->numTimeClusters(); ++tc) {
    // For each cluster, we first store the Ghost cells, then the Copy cells and finally the
    // Interior cells.
    const auto offsetGhost = static_cast<std::size_t>(0);
    const auto offsetCopy = m_ltsTree->layer(LayerIdentifier(HaloType::Ghost, Config(), tc)).size();
    const auto offsetInterior =
        offsetCopy + m_ltsTree->layer(LayerIdentifier(HaloType::Copy, Config(), tc)).size();
    clustersLayerOffset[tc] = LayerOffset{offsetGhost, offsetCopy, offsetInterior};
    clusters[tc + 1] = clusters[tc] + offsetInterior +
                       m_ltsTree->layer(LayerIdentifier(HaloType::Interior, Config(), tc)).size();
  }

  m_meshToClusters.resize(numberOfMeshIds);
  m_meshToLayer.resize(numberOfMeshIds);
  std::size_t cluster = 0;
  std::size_t curClusterElements = 0;
  for (std::size_t cell = 0; cell < numberOfCells; ++cell) {
    if (cell >= clusters[cluster + 1]) {
      curClusterElements = 0;
      ++cluster;
      assert(cluster <= ltsTree->numTimeClusters());
    } else {
      ++curClusterElements;
    }
    if (ltsToMesh[cell] != std::numeric_limits<std::size_t>::max()) {
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
