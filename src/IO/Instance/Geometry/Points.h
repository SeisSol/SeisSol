// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_POINTS_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_POINTS_H_

#include <array>
#include <vector>
namespace seissol::io::instance::geometry {

// equidistant sampling points in GMSH/VTK order

std::vector<std::array<double, 1>> pointsLine(int order);
std::vector<std::array<double, 2>> pointsTriangle(int order);
std::vector<std::array<double, 2>> pointsQuadrangle(int order);
std::vector<std::array<double, 3>> pointsTetrahedron(int order);
std::vector<std::array<double, 3>> pointsHexahedron(int order);

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_POINTS_H_
