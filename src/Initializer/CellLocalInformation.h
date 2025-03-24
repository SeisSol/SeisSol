// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_CELLLOCALINFORMATION_H_
#define SEISSOL_SRC_INITIALIZER_CELLLOCALINFORMATION_H_

#include <Initializer/BasicTypedefs.h>
#include <cstddef>

namespace seissol {

// cell local information
struct CellLocalInformation {
  // types of the faces
  FaceType faceTypes[4];

  // mapping of the neighboring elements to the references element in relation to this element
  int faceRelations[4][2];

  // neighbor config IDs
  unsigned int neighborConfigIds[4];

  // LTS setup
  unsigned short ltsSetup;
};

// cell local information which is not needed during the main iterations, but only during setup and
// (maybe) output
struct SecondaryCellLocalInformation {
  // global mesh ID
  std::size_t globalId;

  // local mesh ID (for interior/copy) or position in the linearized ghost layer
  unsigned int meshId;

  // ids of the face neighbors (in their respective Config LTS tree)
  unsigned int faceNeighborIds[4];

  // ID in layer
  unsigned int layerId;

  // own config ID
  unsigned int configId;

  // unique global id of the time cluster
  unsigned int clusterId;

  // interior, ghost or copy layer cell
  HaloType halo;

  // rank of the own cell
  int rank;

  // rank of all neighboring cells
  int neighborRanks[4];

  // duplicate id of own cell
  char duplicate;
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_CELLLOCALINFORMATION_H_
