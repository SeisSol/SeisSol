// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_

#include <vector>
namespace seissol::io::instance::geometry {

extern const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine4;
extern const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine8;
extern const std::vector<std::vector<std::array<double, 2>>> TriangleRefine4;

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_
