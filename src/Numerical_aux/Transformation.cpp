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
                                                             DenseMatrixView<3, 3> o_Tinv )
{
  real nx = i_normal[0], ny = i_normal[1], nz = i_normal[2];
  real sx = i_tangent1[0], sy = i_tangent1[1], sz = i_tangent1[2];
  real tx = i_tangent2[0], ty = i_tangent2[1], tz = i_tangent2[2];

  o_Tinv(0,0) = nx;
  o_Tinv(1,0) = sx;
  o_Tinv(2,0) = tx;
  o_Tinv(0,1) = ny;
  o_Tinv(1,1) = sy;
  o_Tinv(2,1) = ty;
  o_Tinv(0,2) = nz;
  o_Tinv(1,2) = sz;
  o_Tinv(2,2) = tz;
}

void seissol::transformations::tensor1RotationMatrix( VrtxCoords const i_normal,
                                                      VrtxCoords const i_tangent1,
                                                      VrtxCoords const i_tangent2,
                                                      DenseMatrixView<3, 3> o_T )
{
  real nx = i_normal[0], ny = i_normal[1], nz = i_normal[2];
  real sx = i_tangent1[0], sy = i_tangent1[1], sz = i_tangent1[2];
  real tx = i_tangent2[0], ty = i_tangent2[1], tz = i_tangent2[2];
  
  o_T(0,0) = nx;
  o_T(1,0) = ny;
  o_T(2,0) = nz;
  o_T(0,1) = sx;
  o_T(1,1) = sy;
  o_T(2,1) = sz;
  o_T(0,2) = tx;
  o_T(1,2) = ty;
  o_T(2,2) = tz;
}

void seissol::transformations::inverseSymmetricTensor2RotationMatrix( VrtxCoords const i_normal,
                                                                      VrtxCoords const i_tangent1,
                                                                      VrtxCoords const i_tangent2,
                                                                      DenseMatrixView<6, 6> o_Tinv )
{
  real nx = i_normal[0], ny = i_normal[1], nz = i_normal[2];
  real sx = i_tangent1[0], sy = i_tangent1[1], sz = i_tangent1[2];
  real tx = i_tangent2[0], ty = i_tangent2[1], tz = i_tangent2[2];
  
  o_Tinv(0,0) = nx * nx;
  o_Tinv(1,0) = sx * sx;
  o_Tinv(2,0) = tx * tx;
  o_Tinv(3,0) = nx * sx;
  o_Tinv(4,0) = sx * tx;
  o_Tinv(5,0) = nx * tx;
  o_Tinv(0,1) = ny * ny;
  o_Tinv(1,1) = sy * sy;
  o_Tinv(2,1) = ty * ty;
  o_Tinv(3,1) = ny * sy;
  o_Tinv(4,1) = sy * ty;
  o_Tinv(5,1) = ny * ty;
  o_Tinv(0,2) = nz * nz;
  o_Tinv(1,2) = sz * sz;
  o_Tinv(2,2) = tz * tz;
  o_Tinv(3,2) = nz * sz;
  o_Tinv(4,2) = sz * tz;
  o_Tinv(5,2) = nz * tz;
  o_Tinv(0,3) = 2.0 * ny * nx;
  o_Tinv(1,3) = 2.0 * sy * sx;
  o_Tinv(2,3) = 2.0 * ty * tx;
  o_Tinv(3,3) = ny * sx + nx * sy;
  o_Tinv(4,3) = sy * tx + sx * ty;
  o_Tinv(5,3) = ny * tx + nx * ty;
  o_Tinv(0,4) = 2.0 * nz * ny;
  o_Tinv(1,4) = 2.0 * sz * sy;
  o_Tinv(2,4) = 2.0 * tz * ty;
  o_Tinv(3,4) = nz * sy + ny * sz;
  o_Tinv(4,4) = sz * ty + sy * tz;
  o_Tinv(5,4) = nz * ty + ny * tz;
  o_Tinv(0,5) = 2.0 * nz * nx;
  o_Tinv(1,5) = 2.0 * sz * sx;
  o_Tinv(2,5) = 2.0 * tz * tx;
  o_Tinv(3,5) = nz * sx + nx * sz;
  o_Tinv(4,5) = sz * tx + sx * tz;
  o_Tinv(5,5) = nz * tx + nx * tz;
}

void seissol::transformations::symmetricTensor2RotationMatrix( VrtxCoords const i_normal,
                                                               VrtxCoords const i_tangent1,
                                                               VrtxCoords const i_tangent2,
                                                               DenseMatrixView<6, 6> o_T )
{
  real nx = i_normal[0], ny = i_normal[1], nz = i_normal[2];
  real sx = i_tangent1[0], sy = i_tangent1[1], sz = i_tangent1[2];
  real tx = i_tangent2[0], ty = i_tangent2[1], tz = i_tangent2[2];
  
  o_T(0,0) = nx * nx;
  o_T(1,0) = ny * ny;
  o_T(2,0) = nz * nz;
  o_T(3,0) = ny * nx;
  o_T(4,0) = nz * ny;
  o_T(5,0) = nz * nx;
  o_T(0,1) = sx * sx;
  o_T(1,1) = sy * sy;
  o_T(2,1) = sz * sz;
  o_T(3,1) = sy * sx;
  o_T(4,1) = sz * sy;
  o_T(5,1) = sz * sx;
  o_T(0,2) = tx * tx;
  o_T(1,2) = ty * ty;
  o_T(2,2) = tz * tz;
  o_T(3,2) = ty * tx;
  o_T(4,2) = tz * ty;
  o_T(5,2) = tz * tx;
  o_T(0,3) = 2.0 * nx * sx;
  o_T(1,3) = 2.0 * ny * sy;
  o_T(2,3) = 2.0 * nz * sz;
  o_T(3,3) = ny * sx + nx * sy;
  o_T(4,3) = nz * sy + ny * sz;
  o_T(5,3) = nz * sx + nx * sz;
  o_T(0,4) = 2.0 * sx * tx;
  o_T(1,4) = 2.0 * sy * ty;
  o_T(2,4) = 2.0 * sz * tz;
  o_T(3,4) = sy * tx + sx * ty;
  o_T(4,4) = sz * ty + sy * tz;
  o_T(5,4) = sz * tx + sx * tz;
  o_T(0,5) = 2.0 * nx * tx;
  o_T(1,5) = 2.0 * ny * ty;
  o_T(2,5) = 2.0 * nz * tz;
  o_T(3,5) = ny * tx + nx * ty;
  o_T(4,5) = nz * ty + ny * tz;
  o_T(5,5) = nz * tx + nx * tz;
}



void seissol::transformations::chiTau2XiEtaZeta(unsigned face, double const chiTau[2], double xiEtaZeta[3]) {
  switch (face) {
    case 0:
      xiEtaZeta[0] = chiTau[1];
      xiEtaZeta[1] = chiTau[0];
      xiEtaZeta[2] = 0.0;
      break;
    case 1:
      xiEtaZeta[0] = chiTau[0];
      xiEtaZeta[1] = 0.0;
      xiEtaZeta[2] = chiTau[1];
      break;
    case 2:
      xiEtaZeta[0] = 0.0;
      xiEtaZeta[1] = chiTau[1];
      xiEtaZeta[2] = chiTau[0];
      break;
    case 3:
      xiEtaZeta[0] = 1.0-chiTau[0]-chiTau[1];
      xiEtaZeta[1] = chiTau[0];
      xiEtaZeta[2] = chiTau[1];
      break;
    default:
      break;
  }
}












