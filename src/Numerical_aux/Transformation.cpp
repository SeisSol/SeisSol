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
 * Setup of SeisSol's cell local matrices.
 **/

#include "Transformation.h"
#include <glm/glm.hpp>

void seissol::transformations::tetrahedronReferenceToGlobal( double const v0[3],
                                                             double const v1[3],
                                                             double const v2[3],
                                                             double const v3[3],
                                                             double const xiEtaZeta[3],
                                                             double       xyz[3] ) {
  for (unsigned i = 0; i < 3; ++i) {
    xyz[i] = v0[i] + (v1[i]-v0[i])*xiEtaZeta[0] + (v2[i]-v0[i])*xiEtaZeta[1] + (v3[i]-v0[i])*xiEtaZeta[2];
  }
}

glm::dvec3 seissol::transformations::tetrahedronGlobalToReference( double const       v0[3],
                                                                   double const       v1[3],
                                                                   double const       v2[3],
                                                                   double const      v3[3],
                                                                   glm::dvec3 const& xyz ) {
  // Forward transformation
  glm::dmat4 A( glm::dvec4(v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2], 0.0),
                glm::dvec4(v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2], 0.0),
                glm::dvec4(v3[0]-v0[0], v3[1]-v0[1], v3[2]-v0[2], 0.0),
                glm::dvec4(v0[0], v0[1], v0[2], 1.0)
              );

  double inverseDeterminant = 1.0 / glm::determinant(A);

  glm::dvec4 rhs(xyz[0], xyz[1], xyz[2], 1.0);

  glm::dvec3 xiEtaZeta;
  // Rule of Cramer
  xiEtaZeta[0] = glm::determinant(glm::dmat4(rhs,  A[1], A[2], A[3])) * inverseDeterminant;
  xiEtaZeta[1] = glm::determinant(glm::dmat4(A[0], rhs,  A[2], A[3])) * inverseDeterminant;
  xiEtaZeta[2] = glm::determinant(glm::dmat4(A[0], A[1], rhs,  A[3])) * inverseDeterminant;

  return xiEtaZeta;
}

void seissol::transformations::tetrahedronGlobalToReferenceJacobian( real const i_x[4],
                                                                     real const i_y[4],
                                                                     real const i_z[4],
                                                                     real o_gradXi[3],
                                                                     real o_gradEta[3],
                                                                     real o_gradZeta[3] )
{
  real determinant = i_x[0] * ( i_y[1]*(i_z[3]-i_z[2]) + i_y[2]*(i_z[1]-i_z[3]) + i_y[3]*(i_z[2]-i_z[1]) )
                   + i_x[1] * ( i_y[0]*(i_z[2]-i_z[3]) + i_y[2]*(i_z[3]-i_z[0]) + i_y[3]*(i_z[0]-i_z[2]) )
                   + i_x[2] * ( i_y[0]*(i_z[3]-i_z[1]) + i_y[1]*(i_z[0]-i_z[3]) + i_y[3]*(i_z[1]-i_z[0]) )
                   + i_x[3] * ( i_y[0]*(i_z[1]-i_z[2]) + i_y[1]*(i_z[2]-i_z[0]) + i_y[2]*(i_z[0]-i_z[1]) );                   
  real inverseDeterminant = 1.0 / determinant;
  
  // dxi_dx, dxi_dy, dxi_dz
  o_gradXi[0] = inverseDeterminant * ( i_y[0]*(i_z[2]-i_z[3]) + i_y[2]*(i_z[3]-i_z[0]) + i_y[3]*(i_z[0]-i_z[2]) );
  o_gradXi[1] = inverseDeterminant * ( i_x[0]*(i_z[3]-i_z[2]) + i_x[2]*(i_z[0]-i_z[3]) + i_x[3]*(i_z[2]-i_z[0]) );
  o_gradXi[2] = inverseDeterminant * ( i_x[0]*(i_y[2]-i_y[3]) + i_x[2]*(i_y[3]-i_y[0]) + i_x[3]*(i_y[0]-i_y[2]) );
  
  // deta_dx, deta_dy, deta_dz
  o_gradEta[0] = inverseDeterminant * ( i_y[0]*(i_z[3]-i_z[1]) + i_y[1]*(i_z[0]-i_z[3]) + i_y[3]*(i_z[1]-i_z[0]) );
  o_gradEta[1] = inverseDeterminant * ( i_x[0]*(i_z[1]-i_z[3]) + i_x[1]*(i_z[3]-i_z[0]) + i_x[3]*(i_z[0]-i_z[1]) );
  o_gradEta[2] = inverseDeterminant * ( i_x[0]*(i_y[3]-i_y[1]) + i_x[1]*(i_y[0]-i_y[3]) + i_x[3]*(i_y[1]-i_y[0]) );
  
  // dzeta_dx, dzeta_dx, dzeta_dx
  o_gradZeta[0] = inverseDeterminant * ( i_y[0]*(i_z[1]-i_z[2]) + i_y[1]*(i_z[2]-i_z[0]) + i_y[2]*(i_z[0]-i_z[1]) );
  o_gradZeta[1] = inverseDeterminant * ( i_x[0]*(i_z[2]-i_z[1]) + i_x[1]*(i_z[0]-i_z[2]) + i_x[2]*(i_z[1]-i_z[0]) );
  o_gradZeta[2] = inverseDeterminant * ( i_x[0]*(i_y[1]-i_y[2]) + i_x[1]*(i_y[2]-i_y[0]) + i_x[2]*(i_y[0]-i_y[1]) );
}

void seissol::transformations::inverseTensor1RotationMatrix( VrtxCoords const i_normal,
                                                             VrtxCoords const i_tangent1,
                                                             VrtxCoords const i_tangent2,
                                                             yateto::DenseTensorView<2,real,unsigned>& o_Tinv,
                                                             unsigned row,
                                                             unsigned col )
{
  for (unsigned i = 0; i < 3; ++i) {
    o_Tinv(row+0,col+i) = i_normal[i];
    o_Tinv(row+1,col+i) = i_tangent1[i];
    o_Tinv(row+2,col+i) = i_tangent2[i];
  }
}

void seissol::transformations::tensor1RotationMatrix( VrtxCoords const i_normal,
                                                      VrtxCoords const i_tangent1,
                                                      VrtxCoords const i_tangent2,
                                                      yateto::DenseTensorView<2,real,unsigned>& o_T,
                                                      unsigned row,
                                                      unsigned col )
{
  for (unsigned i = 0; i < 3; ++i) {
    o_T(row+i,col+0) = i_normal[i];
    o_T(row+i,col+1) = i_tangent1[i];
    o_T(row+i,col+2) = i_tangent2[i];
  }
}

void seissol::transformations::inverseSymmetricTensor2RotationMatrix( VrtxCoords const i_normal,
                                                                      VrtxCoords const i_tangent1,
                                                                      VrtxCoords const i_tangent2,
                                                                      yateto::DenseTensorView<2,real,unsigned>& o_Tinv,
                                                                      unsigned row,
                                                                      unsigned col )
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

void seissol::transformations::symmetricTensor2RotationMatrix( VrtxCoords const i_normal,
                                                               VrtxCoords const i_tangent1,
                                                               VrtxCoords const i_tangent2,
                                                               yateto::DenseTensorView<2,real,unsigned>& o_T,
                                                               unsigned row,
                                                               unsigned col )
{  
  real nx = i_normal[0], ny = i_normal[1], nz = i_normal[2];
  real sx = i_tangent1[0], sy = i_tangent1[1], sz = i_tangent1[2];
  real tx = i_tangent2[0], ty = i_tangent2[1], tz = i_tangent2[2];
  
  o_T(row+0,col+0) = nx * nx;
  o_T(row+1,col+0) = ny * ny;
  o_T(row+2,col+0) = nz * nz;
  o_T(row+3,col+0) = ny * nx;
  o_T(row+4,col+0) = nz * ny;
  o_T(row+5,col+0) = nz * nx;
  o_T(row+0,col+1) = sx * sx;
  o_T(row+1,col+1) = sy * sy;
  o_T(row+2,col+1) = sz * sz;
  o_T(row+3,col+1) = sy * sx;
  o_T(row+4,col+1) = sz * sy;
  o_T(row+5,col+1) = sz * sx;
  o_T(row+0,col+2) = tx * tx;
  o_T(row+1,col+2) = ty * ty;
  o_T(row+2,col+2) = tz * tz;
  o_T(row+3,col+2) = ty * tx;
  o_T(row+4,col+2) = tz * ty;
  o_T(row+5,col+2) = tz * tx;
  o_T(row+0,col+3) = 2.0 * nx * sx;
  o_T(row+1,col+3) = 2.0 * ny * sy;
  o_T(row+2,col+3) = 2.0 * nz * sz;
  o_T(row+3,col+3) = ny * sx + nx * sy;
  o_T(row+4,col+3) = nz * sy + ny * sz;
  o_T(row+5,col+3) = nz * sx + nx * sz;
  o_T(row+0,col+4) = 2.0 * sx * tx;
  o_T(row+1,col+4) = 2.0 * sy * ty;
  o_T(row+2,col+4) = 2.0 * sz * tz;
  o_T(row+3,col+4) = sy * tx + sx * ty;
  o_T(row+4,col+4) = sz * ty + sy * tz;
  o_T(row+5,col+4) = sz * tx + sx * tz;
  o_T(row+0,col+5) = 2.0 * nx * tx;
  o_T(row+1,col+5) = 2.0 * ny * ty;
  o_T(row+2,col+5) = 2.0 * nz * tz;
  o_T(row+3,col+5) = ny * tx + nx * ty;
  o_T(row+4,col+5) = nz * ty + ny * tz;
  o_T(row+5,col+5) = nz * tx + nx * tz;
}



void seissol::transformations::chiTau2XiEtaZeta(unsigned face, double const chiTau[2], double xiEtaZeta[3], int sideOrientation) {
  double chiTauTilde[2];

  switch (sideOrientation) {
    case 0:
      chiTauTilde[0] = chiTau[1];
      chiTauTilde[1] = chiTau[0];
      break;
    case 1:
      chiTauTilde[0] = 1.0 - chiTau[0] - chiTau[1];
      chiTauTilde[1] = chiTau[1];
      break;
    case 2:
      chiTauTilde[0] = chiTau[0];
      chiTauTilde[1] = 1.0 - chiTau[0] - chiTau[1];
      break;
    default:
      chiTauTilde[0] = chiTau[0];
      chiTauTilde[1] = chiTau[1];
      break;
  }

  switch (face) {
    case 0:
      xiEtaZeta[0] = chiTauTilde[1];
      xiEtaZeta[1] = chiTauTilde[0];
      xiEtaZeta[2] = 0.0;
      break;
    case 1:
      xiEtaZeta[0] = chiTauTilde[0];
      xiEtaZeta[1] = 0.0;
      xiEtaZeta[2] = chiTauTilde[1];
      break;
    case 2:
      xiEtaZeta[0] = 0.0;
      xiEtaZeta[1] = chiTauTilde[1];
      xiEtaZeta[2] = chiTauTilde[0];
      break;
    case 3:
      xiEtaZeta[0] = 1.0-chiTauTilde[0]-chiTauTilde[1];
      xiEtaZeta[1] = chiTauTilde[0];
      xiEtaZeta[2] = chiTauTilde[1];
      break;
    default:
      break;
  }
}
