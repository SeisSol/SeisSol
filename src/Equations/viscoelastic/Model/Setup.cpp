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

#include <generated_code/sizes.h>
#include <generated_code/init.h>

void getTransposedViscoelasticCoefficientMatrix( real            i_omega,
                                                 unsigned        i_dim,
                                                 unsigned        mech,
                                                 MatrixView      o_M )
{
  unsigned col = 9 + mech * 6;
  switch (i_dim)
  {
    case 0:
      o_M(6, col)     = -i_omega;
      o_M(7, col + 3) = -0.5 * i_omega;
      o_M(8, col + 5) = -0.5 * i_omega;
      break;
      
    case 1:
      o_M(7, col + 1) = -i_omega;
      o_M(6, col + 3) = -0.5 * i_omega;
      o_M(8, col + 4) = -0.5 * i_omega;
      break;
      
    case 2:
      o_M(8, col + 2) = -i_omega;
      o_M(7, col + 4) = -0.5 * i_omega;
      o_M(6, col + 5) = -0.5 * i_omega;
      break;
  }
}

void seissol::model::getTransposedCoefficientMatrix( Material const& i_material,
                                                     unsigned        i_dim,
                                                     real            o_M[seissol::model::AstarT::reals] )
{
  MatrixView M(o_M, seissol::model::AstarT::reals, seissol::model::AstarT::index);  
  // M.setZero();

  seissol::model::getTransposedElasticCoefficientMatrix(i_material, i_dim, M);
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    getTransposedViscoelasticCoefficientMatrix( i_material.omega[mech],
                                                i_dim,
                                                mech,
                                                M );
  }
}

void seissol::model::getTransposedRiemannSolver( seissol::model::Material const&                        local,
                                                 seissol::model::Material const&                        neighbor,
                                                 enum faceType                                          type,
                                                 //real const                                             Atransposed[STAR_NNZ],
                                                 DenseMatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> Flocal,
                                                 DenseMatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> Fneighbor )
{
  real QgodNeighborData[9 * 9];
  real QgodLocalData[9 * 9];
  
  DenseMatrixView<9, 9> QgodNeighbor(QgodNeighborData);
  DenseMatrixView<9, 9> QgodLocal(QgodLocalData);
  
  seissol::model::getTransposedElasticGodunovState(local, neighbor, QgodLocal, QgodNeighbor);
  
  // \todo Generate a kernel for this and use Atransposed instead of the following.
  real tmp[9 * NUMBER_OF_QUANTITIES];
  MatrixView At(tmp, sizeof(tmp)/sizeof(real), &colMjrIndex<9>);
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
    getTransposedViscoelasticCoefficientMatrix(local.omega[mech], 0, mech, At);
    for (unsigned j = 0; j < 6; ++j) {
      unsigned col =  9 + mech * 6 + j;
      for (unsigned i = 0; i < 9; ++i) {
        for (unsigned k = 0; k < 9; ++k) {
          Flocal(i,col) += QgodLocal(i,k) * At(k,col);
          Fneighbor(i,col) += QgodNeighbor(i,k) * At(k,col);
        }
      }
    }
  }
  
  seissol::model::applyBoundaryConditionToElasticFluxSolver(type, Fneighbor.block<9, NUMBER_OF_QUANTITIES>(0, 0));
}

void seissol::model::setMaterial( double* i_materialVal,
                                  int i_numMaterialVals,
                                  seissol::model::Material* o_material )
{
  assert(i_numMaterialVals == 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);
 
  o_material->rho = i_materialVal[0];
  o_material->mu = i_materialVal[1];
  o_material->lambda = i_materialVal[2];
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    o_material->omega[mech] = i_materialVal[3 + 4*mech];
    for (unsigned i = 1; i < 4; ++i) {
      o_material->theta[mech][i-1] = i_materialVal[3 + 4*mech + i];
    }
  }
}


void seissol::model::getFaceRotationMatrix( VrtxCoords const i_normal,
                                            VrtxCoords const i_tangent1,
                                            VrtxCoords const i_tangent2,
                                            DenseMatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> o_T,
                                            DenseMatrixView<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES> o_Tinv )
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

void seissol::model::initializeSpecificLocalData( seissol::model::Material const& material,
                                                  seissol::model::LocalData* localData )
{
  MatrixView sourceMatrix(localData->sourceMatrix, seissol::model::source::reals, seissol::model::source::index);
  sourceMatrix.setZero();

  //       | E_1^T |
  // E^T = |  ...  |
  //       | E_L^T |
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    unsigned offset = 9 + mech * 6;
    real const* theta = material.theta[mech];
    sourceMatrix(offset,     0) = theta[0];
    sourceMatrix(offset + 1, 0) = theta[1];
    sourceMatrix(offset + 2, 0) = theta[1];
    sourceMatrix(offset,     1) = theta[1];
    sourceMatrix(offset + 1, 1) = theta[0];
    sourceMatrix(offset + 2, 1) = theta[1];  
    sourceMatrix(offset,     2) = theta[1];
    sourceMatrix(offset + 1, 2) = theta[1];
    sourceMatrix(offset + 2, 2) = theta[0];  
    sourceMatrix(offset + 3, 3) = theta[2];
    sourceMatrix(offset + 4, 4) = theta[2];
    sourceMatrix(offset + 5, 5) = theta[2];    
  }
  
  // E' = diag(-omega_1 I, ..., -omega_L I)
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    for (unsigned i = 0; i < 6; ++i) {
      unsigned idx = 9 + 6*mech + i;
      sourceMatrix(idx, idx) = -material.omega[mech];
    }
  }
}

void seissol::model::initializeSpecificNeighborData(  seissol::model::Material const&,
                                                      seissol::model::Material const (&)[4],
                                                      seissol::model::NeighborData* )
{
}
