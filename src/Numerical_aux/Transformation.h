/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 **/

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_

#include <yateto.h>
#include <Initializer/typedefs.hpp>
#include <Geometry/MeshDefinition.h>
#include <Eigen/Dense>

namespace seissol {
  namespace transformations {
    /**
     * Calculates the global coordinates from
     * reference tetrahedron coordinates.
     */
    void tetrahedronReferenceToGlobal( double const v0[3],
                                       double const v1[3],
                                       double const v2[3],
                                       double const v3[3],
                                       double const xiEtaZeta[3],
                                       double       xyz[3] );

    Eigen::Vector3d tetrahedronReferenceToGlobal( const Eigen::Vector3d& v0,
                                                  const Eigen::Vector3d& v1,
                                                  const Eigen::Vector3d& v2,
                                                  const Eigen::Vector3d& v3,
                                                  double const xiEtaZeta[3]);

    /**
     * Calculates the reference tetrahedron coordinates from
     * global tetrahedron coordinates.
     */
    Eigen::Vector3d tetrahedronGlobalToReference( double const           v0[3],
                                                  double const           v1[3],
                                                  double const           v2[3],
                                                  double const           v3[3],
                                                  Eigen::Vector3d const& xyz );

    /**
     * Calculates the Jacobian for the coordinate transformation
     * xi(x, y, z), eta(x, y, z), zeta(x, y, z)
     * from a global tetrahedron to the reference tetrahedron.
     **/    
    void tetrahedronGlobalToReferenceJacobian( real const i_x[4],
                                               real const i_y[4],
                                               real const i_z[4],
                                               real o_gradXi[3],
                                               real o_gradEta[3],
                                               real o_gradZeta[3] );

    /**
     * Inverse of Tensor1RotationMatrix().
     **/
    void inverseTensor1RotationMatrix( VrtxCoords const i_normal,
                                       VrtxCoords const i_tangent1,
                                       VrtxCoords const i_tangent2,
                                       yateto::DenseTensorView<2,real,unsigned>& o_Tinv,
                                       unsigned row = 0,
                                       unsigned col = 0 );

    /**
     * Returns a column-major matrix that rotates a first-order tensor
     * (u_x, u_y, u_y) into a new coordinate system aligned with normal
     * and tangents.
     * u' = T*u
     **/
    void tensor1RotationMatrix( VrtxCoords const i_normal,
                                VrtxCoords const i_tangent1,
                                VrtxCoords const i_tangent2,
                                yateto::DenseTensorView<2,real,unsigned>& o_T,
                                unsigned row = 0,
                                unsigned col = 0 );

    /**
     * Inverse of SymmetricTensor2RotationMatrix().
     **/
    template<typename Tmatrix>
    void inverseSymmetricTensor2RotationMatrix( VrtxCoords const i_normal,
                                                VrtxCoords const i_tangent1,
                                                VrtxCoords const i_tangent2,
                                                Tmatrix& o_Tinv,
                                                unsigned row = 0,
                                                unsigned col = 0 )
  {
    real nx = i_normal[0], ny = i_normal[1], nz = i_normal[2];
    real sx = i_tangent1[0], sy = i_tangent1[1], sz = i_tangent1[2];
    real tx = i_tangent2[0], ty = i_tangent2[1], tz = i_tangent2[2];

    o_Tinv(row+0,col+0) = nx * nx;
    o_Tinv(row+1,col+0) = sx * sx;
    o_Tinv(row+2,col+0) = tx * tx;
    o_Tinv(row+3,col+0) = nx * sx;
    o_Tinv(row+4,col+0) = sx * tx;
    o_Tinv(row+5,col+0) = nx * tx;
    o_Tinv(row+0,col+1) = ny * ny;
    o_Tinv(row+1,col+1) = sy * sy;
    o_Tinv(row+2,col+1) = ty * ty;
    o_Tinv(row+3,col+1) = ny * sy;
    o_Tinv(row+4,col+1) = sy * ty;
    o_Tinv(row+5,col+1) = ny * ty;
    o_Tinv(row+0,col+2) = nz * nz;
    o_Tinv(row+1,col+2) = sz * sz;
    o_Tinv(row+2,col+2) = tz * tz;
    o_Tinv(row+3,col+2) = nz * sz;
    o_Tinv(row+4,col+2) = sz * tz;
    o_Tinv(row+5,col+2) = nz * tz;
    o_Tinv(row+0,col+3) = 2.0 * ny * nx;
    o_Tinv(row+1,col+3) = 2.0 * sy * sx;
    o_Tinv(row+2,col+3) = 2.0 * ty * tx;
    o_Tinv(row+3,col+3) = ny * sx + nx * sy;
    o_Tinv(row+4,col+3) = sy * tx + sx * ty;
    o_Tinv(row+5,col+3) = ny * tx + nx * ty;
    o_Tinv(row+0,col+4) = 2.0 * nz * ny;
    o_Tinv(row+1,col+4) = 2.0 * sz * sy;
    o_Tinv(row+2,col+4) = 2.0 * tz * ty;
    o_Tinv(row+3,col+4) = nz * sy + ny * sz;
    o_Tinv(row+4,col+4) = sz * ty + sy * tz;
    o_Tinv(row+5,col+4) = nz * ty + ny * tz;
    o_Tinv(row+0,col+5) = 2.0 * nz * nx;
    o_Tinv(row+1,col+5) = 2.0 * sz * sx;
    o_Tinv(row+2,col+5) = 2.0 * tz * tx;
    o_Tinv(row+3,col+5) = nz * sx + nx * sz;
    o_Tinv(row+4,col+5) = sz * tx + sx * tz;
    o_Tinv(row+5,col+5) = nz * tx + nx * tz;
  }

    /**
     * Returns a column-major matrix that rotates a symmetric second-order
     * tensor, given as vector (u_xx, u_yy, u_zz, u_xy, u_yz, u_xz),
     * into a new coordinate system aligned with normal and tangents.
     * u' = T*u
     **/
    void symmetricTensor2RotationMatrix( VrtxCoords const i_normal,
                                         VrtxCoords const i_tangent1,
                                         VrtxCoords const i_tangent2,
                                         yateto::DenseTensorView<2,real,unsigned>& o_Tinv,
                                         unsigned row = 0,
                                         unsigned col = 0 );

    void chiTau2XiEtaZeta(unsigned face, double const chiTau[2], double xiEtaZeta[3], int sideOrientation = -1);
    void XiEtaZeta2chiTau(unsigned face, double const xiEtaZeta[3], double chiTau[2]);
  }
}

#endif
