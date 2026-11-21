// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_GEOMETRY_FAULTORIENTATION_H_
#define SEISSOL_SRC_CORE_GEOMETRY_FAULTORIENTATION_H_

namespace seissol::initializer::parameters {

/// How the reference point of the faults tells the plus side of a fault face (refPointMethod): as a
/// point on the side of the plus element (Point), or as a direction to compare the normal of the
/// face with (Normal).
enum class RefPointMethod : int { Point = 0, Normal = 1 };

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_CORE_GEOMETRY_FAULTORIENTATION_H_
