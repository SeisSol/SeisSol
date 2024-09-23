// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 */
 
#include "Model/Common.h"
#include <cmath>
#include <iostream>

bool seissol::model::testIfAcoustic(real mu) {
  return std::abs(mu) <= std::numeric_limits<real>::epsilon();
}

void seissol::model::getBondMatrix( VrtxCoords const i_normal,
                                    VrtxCoords const i_tangent1,
                                    VrtxCoords const i_tangent2,
                                    real* o_N )
{
  o_N[0*6 + 0] =   i_normal[0]*i_normal[0]; 
  o_N[0*6 + 1] =   i_normal[1]*i_normal[1];
  o_N[0*6 + 2] =   i_normal[2]*i_normal[2];
  o_N[0*6 + 3] = 2*i_normal[2]*i_normal[1];
  o_N[0*6 + 4] = 2*i_normal[2]*i_normal[0];
  o_N[0*6 + 5] = 2*i_normal[1]*i_normal[0];
  o_N[1*6 + 0] =   i_tangent1[0]*i_tangent1[0]; 
  o_N[1*6 + 1] =   i_tangent1[1]*i_tangent1[1];
  o_N[1*6 + 2] =   i_tangent1[2]*i_tangent1[2];
  o_N[1*6 + 3] = 2*i_tangent1[2]*i_tangent1[1];
  o_N[1*6 + 4] = 2*i_tangent1[2]*i_tangent1[0];
  o_N[1*6 + 5] = 2*i_tangent1[1]*i_tangent1[0];
  o_N[2*6 + 0] =   i_tangent2[0]*i_tangent2[0]; 
  o_N[2*6 + 1] =   i_tangent2[1]*i_tangent2[1];
  o_N[2*6 + 2] =   i_tangent2[2]*i_tangent2[2];
  o_N[2*6 + 3] = 2*i_tangent2[2]*i_tangent2[1];
  o_N[2*6 + 4] = 2*i_tangent2[2]*i_tangent2[0];
  o_N[2*6 + 5] = 2*i_tangent2[1]*i_tangent2[0];
  
  o_N[3*6 + 0] = i_tangent1[0]*i_tangent2[0];
  o_N[3*6 + 1] = i_tangent1[1]*i_tangent2[1];
  o_N[3*6 + 2] = i_tangent1[2]*i_tangent2[2];
  o_N[3*6 + 3] = i_tangent1[1]*i_tangent2[2] + i_tangent1[2]*i_tangent2[1];
  o_N[3*6 + 4] = i_tangent1[0]*i_tangent2[2] + i_tangent1[2]*i_tangent2[0];
  o_N[3*6 + 5] = i_tangent1[1]*i_tangent2[0] + i_tangent1[0]*i_tangent2[1];
  o_N[4*6 + 0] = i_normal[0]*i_tangent2[0];
  o_N[4*6 + 1] = i_normal[1]*i_tangent2[1];
  o_N[4*6 + 2] = i_normal[2]*i_tangent2[2];
  o_N[4*6 + 3] = i_normal[1]*i_tangent2[2] + i_normal[2]*i_tangent2[1];
  o_N[4*6 + 4] = i_normal[0]*i_tangent2[2] + i_normal[2]*i_tangent2[0];
  o_N[4*6 + 5] = i_normal[1]*i_tangent2[0] + i_normal[0]*i_tangent2[1];
  o_N[5*6 + 0] = i_normal[0]*i_tangent1[0];
  o_N[5*6 + 1] = i_normal[1]*i_tangent1[1];
  o_N[5*6 + 2] = i_normal[2]*i_tangent1[2];
  o_N[5*6 + 3] = i_normal[1]*i_tangent1[2] + i_normal[2]*i_tangent1[1];
  o_N[5*6 + 4] = i_normal[0]*i_tangent1[2] + i_normal[2]*i_tangent1[0];
  o_N[5*6 + 5] = i_normal[1]*i_tangent1[0] + i_normal[0]*i_tangent1[1];
}
void seissol::model::getFaceRotationMatrix( Eigen::Vector3d const i_normal,
                                            Eigen::Vector3d const i_tangent1,
                                            Eigen::Vector3d const i_tangent2,
                                            init::T::view::type& o_T,
                                            init::Tinv::view::type& o_Tinv )
{
  VrtxCoords n = {i_normal(0), i_normal(1), i_normal(2)};
  VrtxCoords s = {i_tangent1(0), i_tangent1(1), i_tangent1(2)};
  VrtxCoords t = {i_tangent2(0), i_tangent2(1), i_tangent2(2)};
  getFaceRotationMatrix(n, s, t, o_T, o_Tinv);

}

void seissol::model::getFaceRotationMatrix( VrtxCoords const i_normal,
                                            VrtxCoords const i_tangent1,
                                            VrtxCoords const i_tangent2,
                                            init::T::view::type& o_T,
                                            init::Tinv::view::type& o_Tinv )
{
  o_T.setZero();
  o_Tinv.setZero();
  
  seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, 0, 0);
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, 6, 6);
  
  seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, 0, 0);
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, 6, 6);

#ifdef USE_VISCOELASTIC
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    unsigned const origin = 9 + mech * 6;
    seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, origin, origin);
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, origin, origin);
  }
#elif USE_VISCOELASTIC2
  seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, 9, 9);
#elif USE_POROELASTIC
  //pressure
  o_T(9, 9) = 1;
  o_Tinv(9,9) = 1;
  //fluid velocities
  unsigned origin = 10; 
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, origin, origin);
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, origin, origin);
#endif 
}

