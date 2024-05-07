/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
 * Copyright (c) 2023, Intel corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Point source computation.
 **/

#ifndef SOURCETERM_POINTSOURCE_H_
#define SOURCETERM_POINTSOURCE_H_

#include <Initializer/typedefs.hpp>
#include "SourceTerm/typedefs.hpp"
#include <array>

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
void transformMomentTensor(const real i_localMomentTensor[3][3],
                           const real i_localSolidVelocityComponent[3],
                           real i_localPressureComponent,
                           const real i_localFluidVelocityComponent[3],
                           real strike,
                           real dip,
                           real rake,
                           AlignedArray<real, PointSources::TensorSize>& o_forceComponents);
} // namespace seissol::sourceterm

#endif
