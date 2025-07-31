// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "MeshLayout.h"

#include <Common/Iterator.h>
#include <Geometry/MeshReader.h>
#include <Initializer/BasicTypedefs.h>
#include <Memory/Tree/Colormap.h>
#include <Parallel/MPI.h>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::internal {

std::vector<ClusterMap> layoutCells(const std::vector<std::size_t>& color,
                                    const std::vector<std::size_t>& ghostColor,
                                    const LTSColorMap& colormap,
                                    const geometry::MeshReader& meshReader) {
  // add cells to color clusters
  std::vector<ClusterMap> clusters(colormap.size());

  // layout interior/copy w/o halo
  for (std::size_t i = 0; i < meshReader.getElements().size(); ++i) {
    clusters[color[i]].cellMap.push_back(i);
  }

  // prepare halos
  std::map<std::pair<int, std::size_t>, std::size_t> toLinear;
  std::unordered_map<std::size_t, std::pair<int, std::size_t>> fromLinear;
  std::size_t linearId = 0;
  for (const auto& [rank, cells] : meshReader.getMPINeighbors()) {
    for (const auto [i, cell] : common::enumerate(cells.elements)) {
      clusters[ghostColor[linearId]].cellMap.push_back(linearId);

      toLinear[std::pair<int, std::size_t>{rank, i}] = linearId;
      fromLinear[linearId] = std::pair<int, std::size_t>{rank, i};
      ++linearId;
    }
  }

  // we keep the "old" layout:
  // * ghost and copy colors are matched to each other
  // * copy cells are duplicated, once per ghost color

  // clustered halo
  for (std::size_t i = 0; i < clusters.size(); ++i) {
    if (colormap.argument(i).halo == HaloType::Copy) {
      std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> byColor;
      std::size_t newCount = 0;
      for (std::size_t j = 0; j < clusters[i].cellMap.size(); ++j) {
        const auto cell = clusters[i].cellMap[j];
        std::set<std::pair<int, std::size_t>> neighborColors;
        for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
          // TODO: access linear ghost ID here
          const auto& element = meshReader.getElements()[cell];

          if (element.neighborRanks[f] != MPI::mpi.rank()) {
            const auto ghostLinear = toLinear.at({element.neighborRanks[f], element.mpiIndices[f]});
            const auto colorGhost = ghostColor[ghostLinear];
            neighborColors.emplace(element.neighborRanks[f], colorGhost);
          }
        }
        for (const auto& neighborColor : neighborColors) {
          byColor[neighborColor].push_back(cell);
        }
        newCount += neighborColors.size();
      }

      clusters[i].cellMap.clear();
      clusters[i].cellMap.reserve(newCount);
      for (const auto& [id, cells] : byColor) {
        const auto [rank, color] = id;

        std::vector cellsCopy(cells.begin(), cells.end());
        std::sort(cellsCopy.begin(), cellsCopy.end(), [&](const auto& a, const auto& b) {
          return meshReader.getElements()[a].globalId < meshReader.getElements()[b].globalId;
        });
        clusters[i].cellMap.insert(clusters[i].cellMap.end(), cellsCopy.begin(), cellsCopy.end());
        const auto tag = i + color * colormap.size();
        clusters[i].regions.emplace_back(tag, i, color, cells.size(), rank);
      }
    } else if (colormap.argument(i).halo == HaloType::Ghost) {
      const auto& neighbor = meshReader.getMPINeighbors();
      const auto& ghostlayer = meshReader.getGhostlayerMetadata();
      std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> byColor;
      std::size_t newCount = 0;
      for (std::size_t j = 0; j < clusters[i].cellMap.size(); ++j) {
        const auto cell = clusters[i].cellMap[j];
        const auto dCell = fromLinear[cell];
        const auto& metadata = ghostlayer.at(dCell.first)[dCell.second];
        const auto& localInfo = neighbor.at(dCell.first).elements[dCell.second];
        const auto copyColor = color[localInfo.localElement];
        byColor[{dCell.first, copyColor}].emplace_back(cell);
        newCount += 1;
      }

      clusters[i].cellMap.clear();
      clusters[i].cellMap.reserve(newCount);
      for (const auto& [id, cells] : byColor) {
        const auto [rank, color] = id;

        std::vector cellsCopy(cells.begin(), cells.end());
        std::sort(cellsCopy.begin(), cellsCopy.end(), [&](const auto& a, const auto& b) {
          const auto& dA = fromLinear[a];
          const auto& dB = fromLinear[b];
          const auto aId = ghostlayer.at(dA.first)[dA.second].globalId;
          const auto bId = ghostlayer.at(dB.first)[dB.second].globalId;
          return aId < bId;
        });
        clusters[i].cellMap.insert(clusters[i].cellMap.end(), cellsCopy.begin(), cellsCopy.end());
        const auto tag = i * colormap.size() + color;
        clusters[i].regions.emplace_back(tag, i, color, cells.size(), rank);
      }
    }
  }

  return clusters;
}

std::vector<std::vector<std::size_t>> layoutDR(const std::vector<std::size_t>& color,
                                               const std::vector<std::size_t>& ghostColor,
                                               const LTSColorMap& colormap,
                                               const geometry::MeshReader& meshReader) {
  std::vector<std::vector<std::size_t>> clusters(colormap.size());

  std::map<std::pair<int, std::size_t>, std::size_t> toLinear;
  std::unordered_map<std::size_t, std::pair<int, std::size_t>> fromLinear;
  std::size_t linearId = 0;
  for (const auto& [rank, cells] : meshReader.getMPINeighbors()) {
    for (const auto [i, cell] : common::enumerate(cells.elements)) {
      toLinear[std::pair<int, std::size_t>{rank, i}] = linearId;
      fromLinear[linearId] = std::pair<int, std::size_t>{rank, i};
      ++linearId;
    }
  }

  const auto getColor = [&](std::size_t id, std::size_t idOther, std::size_t sideOther) {
    if (id == std::numeric_limits<std::size_t>::max()) {
      const auto rank = meshReader.getElements()[idOther].neighborRanks[sideOther];
      const auto index = meshReader.getElements()[idOther].mpiIndices[sideOther];

      return ghostColor[toLinear.at(std::pair<int, std::size_t>{rank, index})];
    } else {
      return color[id];
    }
  };

  const auto colorAdjust = [&](std::size_t color, HaloType halo) {
    auto id = colormap.argument(color);
    id.halo = halo;
    return colormap.colorId(id);
  };

  const auto colorCompare = [&](std::size_t a, std::size_t b) {
    return colorAdjust(a, HaloType::Interior) == colorAdjust(b, HaloType::Interior);
  };

  const auto determineColor = [&](std::size_t a, std::size_t b) {
    auto idA = colormap.argument(a);
    auto idB = colormap.argument(b);

    if (idA.halo == HaloType::Ghost || idB.halo == HaloType::Ghost) {
      return colorAdjust(a, HaloType::Copy);
    } else {
      return colorAdjust(a, HaloType::Ghost);
    }
  };

  for (const auto [i, fault] : common::enumerate(meshReader.getFault())) {
    const auto plusColor = getColor(fault.element, fault.neighborElement, fault.neighborSide);
    const auto minusColor = getColor(fault.neighborElement, fault.element, fault.neighborSide);

    if (!colorCompare(plusColor, minusColor)) {
      logError() << "";
    }

    clusters[determineColor(plusColor, minusColor)].emplace_back(i);
  }

  return clusters;
}

} // namespace seissol::initializer::internal
