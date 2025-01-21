// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_SOURCETERM_POINTSOURCE_H_
#define SEISSOL_SRC_SOURCETERM_POINTSOURCE_H_

#include "Initializer/Typedefs.h"
#include "SourceTerm/Typedefs.h"

namespace seissol::sourceterm {
/** The local moment tensor shall be transformed into the global coordinate system.
 *
 * The second order tensor (matrix) can be understood as a transform
 * on a vector, e.g. p_L = T_L * q_L. (Let L = Local moment tensor, M = Moment tensor.)
 * We are looking for the transformed tensor p_M = T_M * q_M, i.e.
 * the local moment tensor rotated by strike, dip, and rake.
 * Assume x_L = R * x_M, where R is an orthogonal matrix. Then
 * p_M = R^T * p_L = R^T * T_L * R * q_M and hence
 *   T_M = R^T * T_L * R.
 * Thus, the rotation matrix R is the transformation from the global (x,y,z)
 * coordinate system to local (fault plane) coordinate system and is obtained
 * by the successive rotations  strike (s) -> dip (d) -> rake (l).
 *
 *                   |  cos l  sin l    | | 1               |  |  cos s -sin s    |
 * R_l * R_d * R_s = | -sin l  cos l    | |    cos d -sin d |  |  sin s  cos s    |
 *                   |                1 | |    sin d  cos d |  |                1 |
 *
 **/
void transformMomentTensor(
    const real localMomentTensor[3][3],
    const real localSolidVelocityComponent[3],
    real localPressureComponent,
    const real localFluidVelocityComponent[3],
    real strike,
    real dip,
    real rake,
    seissol::memory::AlignedArray<real, PointSources::TensorSize>& forceComponents);
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_POINTSOURCE_H_
