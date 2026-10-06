// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FAULTGEOMETRY_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FAULTGEOMETRY_H_

#include "Geometry/CellGeometry.h"
#include "Geometry/FaceTransform.h"

#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::geometry {
class MeshReader;
} // namespace seissol::geometry

namespace seissol::dr {

/// The quadrature point `point` of a fault face, in the coordinates of its plus face.
auto quadraturePoint(std::size_t point) -> geometry::FaceTransform::FaceVectorT;

/**
 * The quadrature point a point of the arrays of a fault face stands for, which interleave the
 * fused simulations and are padded to the vector width; none for the padding.
 */
auto quadraturePointOf(std::size_t paddedPoint) -> std::optional<std::size_t>;

/**
 * The frame of the fault face `faultId` at each of its quadrature points (see
 * geometry::faultFrameAt). On a plane face it is the frame of the fault at every point.
 */
auto faultFramesAtPoints(std::size_t faultId, const geometry::MeshReader& mesh)
    -> std::vector<geometry::FaultFrame>;

} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FAULTGEOMETRY_H_
