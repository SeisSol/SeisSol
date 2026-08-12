// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "Initializer/Clustering/ClusterSmoother.h"

#include "Common/Constants.h"
#include "Geometry/PUMLReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Parallel/MPI.h"

#include <PUML/Downward.h>
#include <PUML/Upward.h>
#include <algorithm>
#include <cstddef>
#include <mpi.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::initializer {

FaceType decodeFaceType(const void* boundaryCond,
                        std::size_t cell,
                        std::uint8_t face,
                        parameters::BoundaryFormat boundaryFormat,
                        const FaceMap& faceMap) {
  const auto tag = seissol::geometry::decodeBoundary(boundaryCond, cell, face, boundaryFormat);
  const auto faceType = faceMap.at(tag);
  if (!faceType.has_value()) {
    logError() << "Read invalid face tag during cluster determining:" << tag;
  }
  return faceType.value();
}

// PUML's traversal API works in `unsigned int` face ids and `int` local cell ids, and MPI counts
// are `int`. Those types appear verbatim wherever this file talks to either library; everything
// in between uses `std::size_t`.
ClusterSmoother::ClusterSmoother(const geometry::PumlMesh& mesh,
                                 parameters::BoundaryFormat boundaryFormat,
                                 const FaceMap& faceMap)
    : mesh_(&mesh), boundaryFormat_(boundaryFormat), faceMap_(&faceMap) {
  const auto& cells = mesh_->cells();
  const auto& faces = mesh_->faces();
  const void* boundaryCond = mesh_->cellData(1);

  std::unordered_map<int, std::vector<std::size_t>> rankToSharedFacesPre;
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    unsigned int faceids[Cell::NumFaces]{};
    bool atBoundary = false;
    PUML::Downward::faces(*mesh_, cells[cell], faceids);
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto boundary = decodeFaceType(boundaryCond, cell, f, boundaryFormat_, *faceMap_);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbors later
        const auto& face = faces.at(faceids[f]);
        if (face.isShared()) {
          rankToSharedFacesPre[face.shared()[0]].push_back(faceids[f]);
          localFaceIdToLocalCellId_[faceids[f]] = cell;
          atBoundary = true;
        }
      }
    }
    if (atBoundary) {
      boundaryCells_.emplace_back(cell);
    }
  }

  for (auto& sharedFaces : rankToSharedFacesPre) {
    std::sort(sharedFaces.second.begin(),
              sharedFaces.second.end(),
              [&](unsigned int a, unsigned int b) { return faces[a].gid() < faces[b].gid(); });
  }

  rankToSharedFaces_ =
      decltype(rankToSharedFaces_)(rankToSharedFacesPre.begin(), rankToSharedFacesPre.end());

  for (std::size_t ex = 0; ex < rankToSharedFaces_.size(); ++ex) {
    const auto& exchange = rankToSharedFaces_[ex];
    for (std::size_t i = 0; i < exchange.second.size(); ++i) {
      sharedFaceToExchangeId_[exchange.second[i]] = {ex, i};
    }
  }

  // allocate the exchange buffers once; every sweep reuses them
  const auto numExchanges = rankToSharedFaces_.size();
  requests_.resize(2 * numExchanges);
  ghost_.resize(numExchanges);
  copy_.resize(numExchanges);
  for (std::size_t ex = 0; ex < numExchanges; ++ex) {
    const auto exchangeSize = rankToSharedFaces_[ex].second.size();
    ghost_[ex].resize(exchangeSize);
    copy_[ex].resize(exchangeSize);
  }
}

std::size_t ClusterSmoother::relaxOnce(std::vector<std::size_t>& clusterIds,
                                       const SmoothingRule& rule) {
  std::size_t numberOfReductions = 0;

  const auto& cells = mesh_->cells();
  const auto& faces = mesh_->faces();
  const void* boundaryCond = mesh_->cellData(1);

#pragma omp parallel for reduction(+ : numberOfReductions)
  for (std::size_t cell = 0; cell < cells.size(); ++cell) {
    std::size_t timeCluster = clusterIds[cell];

    unsigned int faceids[Cell::NumFaces]{};
    PUML::Downward::faces(*mesh_, cells[cell], faceids);
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto boundary = decodeFaceType(boundaryCond, cell, f, boundaryFormat_, *faceMap_);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbors later
        const auto& face = faces.at(faceids[f]);
        if (!face.isShared()) {
          int cellIds[2];
          PUML::Upward::cells(*mesh_, face, cellIds);

          // PUML hands back signed local ids; widening to the index type is explicit here
          const int neighborCell = (cellIds[0] == static_cast<int>(cell)) ? cellIds[1] : cellIds[0];
          const auto otherTimeCluster = clusterIds[static_cast<std::size_t>(neighborCell)];

          // the rule yields a non-negative index difference, so the sum cannot wrap
          const auto admissible =
              otherTimeCluster + static_cast<std::size_t>(rule.differenceFor(boundary));

          if (timeCluster > admissible) {
            timeCluster = admissible;
            ++numberOfReductions;
          }
        }
      }
    }
    clusterIds[cell] = timeCluster;
  }

  const auto numExchanges = rankToSharedFaces_.size();

  for (std::size_t ex = 0; ex < numExchanges; ++ex) {
    const auto& exchange = rankToSharedFaces_[ex];
    const auto exchangeSize = exchange.second.size();

    for (std::size_t n = 0; n < exchangeSize; ++n) {
      copy_[ex][n] = clusterIds[localFaceIdToLocalCellId_[exchange.second[n]]];
    }
    MPI_Isend(copy_[ex].data(),
              exchangeSize,
              Mpi::castToMpiType<std::size_t>(),
              exchange.first,
              0,
              seissol::Mpi::mpi.comm(),
              &requests_[ex]);
    MPI_Irecv(ghost_[ex].data(),
              exchangeSize,
              Mpi::castToMpiType<std::size_t>(),
              exchange.first,
              0,
              seissol::Mpi::mpi.comm(),
              &requests_[numExchanges + ex]);
  }

  MPI_Waitall(2 * numExchanges, requests_.data(), MPI_STATUSES_IGNORE);

#pragma omp parallel for reduction(+ : numberOfReductions)
  for (std::size_t bcell = 0; bcell < boundaryCells_.size(); ++bcell) {
    const auto cell = boundaryCells_[bcell];
    std::size_t& timeCluster = clusterIds[cell];

    unsigned int faceids[Cell::NumFaces]{};
    PUML::Downward::faces(*mesh_, cells[cell], faceids);
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto boundary = decodeFaceType(boundaryCond, cell, f, boundaryFormat_, *faceMap_);
      // Continue for regular, dynamic rupture, and periodic boundary cells
      if (isInternalFaceType(boundary)) {
        // We treat MPI neighbors later
        const auto& face = faces.at(faceids[f]);
        if (face.isShared()) {
          const auto pos = sharedFaceToExchangeId_.at(faceids[f]);
          const auto otherTimeCluster = ghost_[pos.first][pos.second];

          const auto admissible =
              otherTimeCluster + static_cast<std::size_t>(rule.differenceFor(boundary));

          if (timeCluster > admissible) {
            timeCluster = admissible;
            ++numberOfReductions;
          }
        }
      }
    }
  }

  return numberOfReductions;
}

std::size_t ClusterSmoother::relax(std::vector<std::size_t>& clusterIds,
                                   const SmoothingRule& rule,
                                   MPI_Comm comm) {
  std::size_t totalNumberOfReductions = 0;
  std::size_t globalNumberOfReductions = 0;
  do {
    std::size_t localNumberOfReductions = relaxOnce(clusterIds, rule);

    MPI_Allreduce(&localNumberOfReductions,
                  &globalNumberOfReductions,
                  1,
                  Mpi::castToMpiType<std::size_t>(),
                  MPI_SUM,
                  comm);
    totalNumberOfReductions += globalNumberOfReductions;
  } while (globalNumberOfReductions > 0);
  return totalNumberOfReductions;
}

} // namespace seissol::initializer
