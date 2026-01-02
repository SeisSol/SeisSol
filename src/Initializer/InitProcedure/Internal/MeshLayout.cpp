// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "MeshLayout.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Geometry/MeshReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Tree/Colormap.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <utility>
#include <utils/logger.h>
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

  for (std::size_t i = 0; i < meshReader.linearGhostlayer().size(); ++i) {
    clusters[ghostColor[i]].cellMap.push_back(i);
  }

  // we keep the "old" layout:
  // * ghost and copy colors are matched to each other
  // * copy cells are duplicated, once per ghost color

  const auto colorAdjust = [&](std::size_t color, HaloType halo) {
    auto id = colormap.argument(color);
    id.halo = halo;
    return colormap.colorId(id);
  };

  // clustered halo
  for (std::size_t i = 0; i < clusters.size(); ++i) {
    if (colormap.argument(i).halo == HaloType::Copy) {
      std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> byColor;
      std::size_t newCount = 0;
      for (std::size_t j = 0; j < clusters[i].cellMap.size(); ++j) {
        const auto cell = clusters[i].cellMap[j];
        std::set<std::pair<int, std::size_t>> neighborColors;
        for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
          const auto& element = meshReader.getElements()[cell];

          if (element.neighborRanks[f] != seissol::Mpi::mpi.rank()) {
            const auto ghostLinear = meshReader.toLinearGhostlayer().at(
                {element.neighborRanks[f], element.mpiIndices[f]});
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
      clusters[i].regions.reserve(byColor.size());
      for (const auto& [id, cells] : byColor) {
        const auto [rank, color] = id;

        std::vector cellsCopy(cells.begin(), cells.end());
        std::sort(cellsCopy.begin(), cellsCopy.end(), [&](const auto& a, const auto& b) {
          return meshReader.getElements()[a].globalId < meshReader.getElements()[b].globalId;
        });
        clusters[i].cellMap.insert(clusters[i].cellMap.end(), cellsCopy.begin(), cellsCopy.end());
        const auto tag = colorAdjust(i, HaloType::Interior) +
                         colorAdjust(color, HaloType::Interior) * colormap.size();
        clusters[i].regions.emplace_back(tag, i, color, cells.size(), rank);
      }
    } else if (colormap.argument(i).halo == HaloType::Ghost) {
      const auto& neighbor = meshReader.getMPINeighbors();
      const auto& ghostlayer = meshReader.getGhostlayerMetadata();

      const auto globalId = [&](std::size_t id) {
        // take the first index here; linearGhostLayer will not contain empty vectors
        const auto& dId = meshReader.linearGhostlayer()[id];
        return ghostlayer.at(dId.rank)[dId.inRankIndices[0]].globalId;
      };

      std::map<std::pair<int, std::size_t>, std::vector<std::size_t>> byColor;
      std::size_t newCount = 0;
      for (std::size_t j = 0; j < clusters[i].cellMap.size(); ++j) {
        const auto cell = clusters[i].cellMap[j];
        const auto& linearCell = meshReader.linearGhostlayer()[cell];

        std::set<std::pair<int, std::size_t>> neighborColors;
        for (const auto index : linearCell.inRankIndices) {
          const auto& localInfo = neighbor.at(linearCell.rank).elements[index];
          const auto copyColor = color[localInfo.localElement];

          neighborColors.emplace(linearCell.rank, copyColor);
        }
        for (const auto& neighborColor : neighborColors) {
          byColor[neighborColor].push_back(cell);
        }
        newCount += neighborColors.size();
      }

      clusters[i].cellMap.clear();
      clusters[i].cellMap.reserve(newCount);
      clusters[i].regions.reserve(byColor.size());
      for (const auto& [id, cells] : byColor) {
        const auto [rank, color] = id;

        std::vector cellsCopy(cells.begin(), cells.end());
        std::sort(cellsCopy.begin(), cellsCopy.end(), [&](const auto& a, const auto& b) {
          return globalId(a) < globalId(b);
        });
        clusters[i].cellMap.insert(clusters[i].cellMap.end(), cellsCopy.begin(), cellsCopy.end());
        const auto tag = colorAdjust(i, HaloType::Interior) +
                         colorAdjust(color, HaloType::Interior) * colormap.size();
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

  const auto getColor = [&](const auto& id, const auto& idOther, std::size_t sideOther) {
    if (id.hasValue()) {
      return color[id.value()];
    } else {
      const auto rank = meshReader.getElements()[idOther.value()].neighborRanks[sideOther];
      const auto index = meshReader.getElements()[idOther.value()].mpiIndices[sideOther];

      return ghostColor[meshReader.toLinearGhostlayer().at(
          std::pair<int, std::size_t>{rank, index})];
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
      return colorAdjust(a, HaloType::Interior);
    }
  };

  for (const auto [i, fault] : common::enumerate(meshReader.getFault())) {
    const auto plusColor = getColor(fault.element, fault.neighborElement, fault.neighborSide);
    const auto minusColor = getColor(fault.neighborElement, fault.element, fault.side);

    if (!colorCompare(plusColor, minusColor)) {
      logError()
          << "The element setup differs around a fault face (currently, SeisSol requires both "
             "sides of a fault face to use the same configuration and LTS cluster). Aborting.";
    }

    clusters[determineColor(plusColor, minusColor)].emplace_back(i);
  }

  return clusters;
}

} // namespace seissol::initializer::internal
