// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_CELLGEOMETRY_H_
#define SEISSOL_SRC_GEOMETRY_CELLGEOMETRY_H_

#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"

#include <cstddef>
#include <memory>

namespace seissol::geometry {

class MeshReader;

/**
 * The transform of a mesh cell, whatever shape the mesh gives it: the affine one where the mesh
 * is straight-sided, the isoparametric one through its nodes where it carries curved cells (see
 * MeshReader::geometryOrder). Everything that maps between the reference cell and space for a
 * mesh cell should go through here, so that a curved mesh reaches it without further ado.
 */
auto cellTransformOf(std::size_t id, const MeshReader& mesh) -> std::unique_ptr<CellTransform>;

/// the face `side` of a mesh cell, as seen from that cell (orientation Local) or from its
/// neighbour (see FaceOrientation)
auto faceTransformOf(std::size_t id,
                     std::size_t side,
                     const MeshReader& mesh,
                     FaceOrientation orientation = FaceOrientation::Local)
    -> std::unique_ptr<FaceTransform>;

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_CELLGEOMETRY_H_
