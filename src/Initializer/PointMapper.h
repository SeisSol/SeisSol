// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_POINTMAPPER_H_
#define SEISSOL_SRC_INITIALIZER_POINTMAPPER_H_

#include "Geometry/MeshReader.h"
#include <Eigen/Dense>

namespace seissol::initializer {
/** Finds the tetrahedrons that contain the points.
 *  In "contained" we save if the point source is contained in the mesh.
 *  We use short here as bool. For MPI use cleanDoubles afterwards.
 */
void findMeshIds(const Eigen::Vector3d* points,
                 const seissol::geometry::MeshReader& mesh,
                 std::size_t numPoints,
                 short* contained,
                 std::size_t* meshId,
                 double tolerance = 0);

void findMeshIds(const Eigen::Vector3d* points,
                 const std::vector<Vertex>& vertices,
                 const std::vector<Element>& elements,
                 std::size_t numPoints,
                 short* contained,
                 std::size_t* meshIds,
                 double tolerance = 0);

void cleanDoubles(short* contained, std::size_t numPoints);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_POINTMAPPER_H_
