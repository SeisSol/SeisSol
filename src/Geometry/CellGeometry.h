// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_CELLGEOMETRY_H_
#define SEISSOL_SRC_GEOMETRY_CELLGEOMETRY_H_

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshDefinition.h"

#include <array>
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

/**
 * The face of the fault `faultId` in the parametrisation of its plus side, from whichever of the
 * two cells sharing it is on this rank (as AffineFaceTransform::fromMeshFault does), curved where
 * the mesh carries curved cells.
 */
auto faultFaceTransformOf(std::size_t faultId, const MeshReader& mesh)
    -> std::unique_ptr<FaceTransform>;

/// The frame of a fault at a point of it
struct FaultFrame {
  /// pointing from the plus side to the minus side
  CellTransform::VectorEigenT normal;
  CellTransform::VectorEigenT tangent1;
  CellTransform::VectorEigenT tangent2;
  /// the surface Jacobian there, FaceTransform::surfaceJacobian
  double surfaceJacobian{};
};

/**
 * The frame of a fault at the point `point` of its face `face` (see faultFaceTransformOf). A plane
 * face has the frame of the fault (Fault::normal, Fault::tangent1, Fault::tangent2) at every point.
 * On a curved face the normal turns, and the first tangent is the one of the fault turned into the
 * plane orthogonal to the normal there -- which the two cells sharing the face, possibly on two
 * ranks, both come to.
 */
auto faultFrameAt(const FaceTransform& face,
                  const Fault& fault,
                  const FaceTransform::FaceVectorT& point) -> FaultFrame;

/**
 * How thin a cell is where it is thinnest, relative to the straight-sided cell through its
 * vertices: the smallest singular value of its Jacobian over the one of that straight-sided cell,
 * the smallest over the nodes of a lattice in the cell. It is one for a straight-sided cell, and
 * what the time step of a curved cell is scaled by: where a cell is squeezed, its metric -- and
 * with it the speed a wave crosses the reference cell at -- grows by the inverse.
 */
auto relativeThickness(const CellTransform& transform,
                       const std::array<CellTransform::VectorEigenT, Cell::NumVertices>& vertices)
    -> double;

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_CELLGEOMETRY_H_
