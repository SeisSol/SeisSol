// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_

#include "Common/Constants.h"
#include "Geometry/MeshReader.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Descriptor/LTS.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace seissol::initializer::internal {

/**
 * Derives the storage requirements of a single cell. Exposed for testing; see the definition for
 * the encoding of the result.
 */
LtsSetup getLtsSetup(const CellLocalInformation& ownPrimary,
                     const SecondaryCellLocalInformation& ownSecondary,
                     const std::array<std::uint64_t, Cell::NumFaces>& neighborClusters);

void deriveLtsSetups(const MeshLayout& layout, LTS::Storage& storage);

/**
 * In which parametrisation the traces of a face are written that a cell shares with a neighbor
 * (CellLocalInformation::traceSides), each given by its global ID and its number of the face: 0,
 * the cell's own, where the cell comes first, and 1, the neighbor's, otherwise. The two cells of a
 * face come to opposite answers, also where a cell is its own neighbor across two of its faces.
 */
constexpr std::uint8_t traceSideOf(std::size_t globalId,
                                   std::size_t face,
                                   std::size_t neighborGlobalId,
                                   std::size_t neighborFace) {
  return std::make_pair(globalId, face) > std::make_pair(neighborGlobalId, neighborFace) ? 1 : 0;
}

/**
 * Decides, for every face of every cell of this rank, which parametrisation of the face the face
 * traces are written in (CellLocalInformation::traceSides): the plus side of a fault keeps its own
 * and the minus side takes it, and of the two cells sharing any other face the one with the
 * smaller global ID keeps its own (traceSideOf). Needs the face neighbors and their global IDs, the
 * ghost layer included.
 */
void deriveTraceSides(const geometry::MeshReader& meshReader, LTS::Storage& storage);
} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_
