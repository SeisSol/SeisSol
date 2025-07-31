// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_CELLLOCALINFORMATION_H_
#define SEISSOL_SRC_INITIALIZER_CELLLOCALINFORMATION_H_

#include <Common/Constants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/LtsSetup.h>
#include <Memory/Tree/Backmap.h>
#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol {

// cell local information
struct CellLocalInformation {
  // types of the faces
  std::array<FaceType, Cell::NumFaces> faceTypes;

  // mapping of the neighboring elements to the references element in relation to this element
  std::array<std::array<uint8_t, 2>, Cell::NumFaces> faceRelations;

  // neighbor config IDs
  std::array<std::uint32_t, 4> neighborConfigIds;

  // LTS setup
  LtsSetup ltsSetup;
};

// cell local information which is not needed during the main iterations, but only during setup and
// (maybe) output
struct SecondaryCellLocalInformation {
  // global mesh ID
  std::size_t globalId;

  // local mesh ID (for interior/copy) or position in the linearized ghost layer
  std::size_t meshId;

  // global ids of the face neighbors
  std::array<std::size_t, 4> faceNeighborIds;

  // full storage positions of the face neighbors
  std::array<initializer::StoragePosition, 4> faceNeighborPositions;

  // ID in layer
  std::size_t layerId;

  // own config ID
  std::uint32_t configId;

  // unique global id of the time cluster
  std::uint64_t clusterId;

  // interior, ghost or copy layer cell
  HaloType halo;

  // rank of the own cell
  int rank;

  // rank of all neighboring cells
  std::array<int, 4> neighborRanks;

  // duplicate id of own cell
  uint8_t duplicate;
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_CELLLOCALINFORMATION_H_
