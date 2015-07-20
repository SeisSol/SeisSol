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

#include <Model/Setup.h>

#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>


void getTransposedViscoelasticCoefficientMatrix( real            i_omega,
                                                 unsigned        i_dim,
                                                 MatrixView<9,6> o_M )
{
  o_M.setZero();
  
  switch (i_dim)
  {
    case 0:
      o_M(6, 0) = -i_omega;
      o_M(7, 3) = -0.5 * i_omega;
      o_M(8, 5) = -0.5 * i_omega;
      break;
      
    case 1:
      o_M(7, 1) = -i_omega;
      o_M(6, 3) = -0.5 * i_omega;
      o_M(8, 4) = -0.5 * i_omega;
      break;
      
    case 2:
      o_M(8, 2) = -i_omega;
      o_M(7, 4) = -0.5 * i_omega;
      o_M(6, 5) = -0.5 * i_omega;
      break;
  }
}

void seissol::model::getTransposedCoefficientMatrix( Material const& i_material,
                                                     unsigned        i_dim,
                                                     real            o_M[STAR_NNZ] )
{
  assert(STAR_NNZ == NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES);
  MatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> M(o_M);
  
  M.setZero();
  seissol::model::getTransposedElasticCoefficientMatrix(i_material, i_dim, M.block<9,9>(0,0));
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    getTransposedViscoelasticCoefficientMatrix( i_material.omega[mech],
                                                i_dim,
                                                M.block<9,6>(0, 9 + mech * 6) );
  }
}

void seissol::model::getTransposedRiemannSolver( seissol::model::Material const&                        local,
                                                 seissol::model::Material const&                        neighbor,
                                                 enum faceType                                          type,
                                                 //real const                                             Atransposed[STAR_NNZ],
                                                 MatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> Flocal,
                                                 MatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> Fneighbor )
{
  real QgodNeighborData[9 * 9];
  real QgodLocalData[9 * 9];
  
  MatrixView<9, 9> QgodNeighbor(QgodNeighborData);
  MatrixView<9, 9> QgodLocal(QgodLocalData);
  
  seissol::model::getTransposedElasticGodunovState(local, neighbor, QgodLocal, QgodNeighbor);
  
  // \todo Generate a kernel for this and use Atransposed instead of the following.
  real tmp[9 * 9];
  MatrixView<9, 9> At(tmp);
  seissol::model::getTransposedElasticCoefficientMatrix(local, 0, At);
  Flocal.setZero();
  Fneighbor.setZero();
  for (unsigned j = 0; j < 9; ++j) {
    for (unsigned i = 0; i < 9; ++i) {
      for (unsigned k = 0; k < 9; ++k) {
        Flocal(i,j) += QgodLocal(i,k) * At(k,j);
        Fneighbor(i,j) += QgodNeighbor(i,k) * At(k,j);
      }
    }
  }
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    unsigned origin = 9 + mech * 6;
    MatrixView<9, 6> At_l(tmp);
    getTransposedViscoelasticCoefficientMatrix(local.omega[mech], 0, At_l);
    for (unsigned j = 0; j < 6; ++j) {
      for (unsigned i = 0; i < 9; ++i) {
        for (unsigned k = 0; k < 9; ++k) {
          Flocal(i,origin + j) += QgodLocal(i,k) * At_l(k,j);
          Fneighbor(i,origin + j) += QgodNeighbor(i,k) * At_l(k,j);
        }
      }
    }
  }
  
  seissol::model::applyBoundaryConditionToElasticFluxSolver(type, Fneighbor.block<NUMBER_OF_QUANTITIES, 9>(0, 0));
}

void seissol::model::setMaterial( double* i_materialVal,
                                  int i_numMaterialVals,
                                  seissol::model::Material* o_material )
{
  assert(i_numMaterialVals == 3);
 
  o_material->rho = i_materialVal[0];
  o_material->mu = i_materialVal[1];
  o_material->lambda = i_materialVal[2];  
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    o_material->omega[mech] = i_materialVal[4 + 4*mech];
  }
}


void seissol::model::getFaceRotationMatrix( VrtxCoords const i_normal,
                                            VrtxCoords const i_tangent1,
                                            VrtxCoords const i_tangent2,
                                            MatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> o_T,
                                            MatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> o_Tinv )
{
  o_T.setZero();
  o_Tinv.setZero();
  
  seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T.block<6,6>(0, 0));
  seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv.block<6,6>(0, 0));
  
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T.block<3,3>(6, 6));
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv.block<3,3>(6, 6));
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    unsigned origin = 9 + mech * 6;
    seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T.block<6,6>(origin, origin));
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv.block<6,6>(origin, origin));
  }
}
