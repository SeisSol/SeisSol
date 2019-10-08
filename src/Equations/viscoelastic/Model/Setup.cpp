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

#include <yateto/TensorView.h>
#include <generated_code/init.h>

template<typename T>
void getTransposedViscoelasticCoefficientMatrix( real            i_omega,
                                                 unsigned        i_dim,
                                                 unsigned        mech,
                                                 T&              o_M )
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

template<typename T>
void getTransposedSourceCoefficientTensor( seissol::model::Material const& material, T& sourceMatrix )
{
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

void seissol::model::getTransposedCoefficientMatrix( Material const& i_material,
                                                     unsigned        i_dim,
                                                     init::star::view<0>::type& AT )
{
  seissol::model::getTransposedElasticCoefficientMatrix(i_material, i_dim, AT);

  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    getTransposedViscoelasticCoefficientMatrix( i_material.omega[mech],
                                                i_dim,
                                                mech,
                                                AT );
  }
}

void seissol::model::getPlaneWaveOperator(  Material const& material,
                                            double const n[3],
                                            std::complex<real> Mdata[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES] )
{
  yateto::DenseTensorView<2,std::complex<real>> M(Mdata, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});
  M.setZero();

  real data[NUMBER_OF_QUANTITIES * NUMBER_OF_QUANTITIES];
  yateto::DenseTensorView<2,real> Coeff(data, {NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES});

  for (unsigned d = 0; d < 3; ++d) {
    Coeff.setZero();
    getTransposedCoefficientMatrix(material, d, Coeff);

    for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
        M(i,j) += n[d] * Coeff(j,i);
      }
    }
  }

  getTransposedSourceCoefficientTensor(material, Coeff);

  for (unsigned i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
    for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; ++j) {
      M(i,j) -= std::complex<real>(0.0, Coeff(j,i));
    }
  }
}

void seissol::model::getTransposedGodunovState( Material const&                   local,
                                                Material const&                   neighbor,
                                                enum ::faceType                   faceType,
                                                init::QgodLocal::view::type&      QgodLocal,
                                                init::QgodNeighbor::view::type&   QgodNeighbor )
{
  seissol::model::getTransposedElasticGodunovState(local, neighbor, faceType, QgodLocal, QgodNeighbor);
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
                                            init::T::view::type& o_T,
                                            init::Tinv::view::type& o_Tinv )
{
  o_T.setZero();
  o_Tinv.setZero();
  
  seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T);
  seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv);
  
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, 6, 6);
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, 6, 6);
  
  for (unsigned mech = 0; mech < NUMBER_OF_RELAXATION_MECHANISMS; ++mech) {
    unsigned origin = 9 + mech * 6;
    seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_T, origin, origin);
    seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, o_Tinv, origin, origin);
  }
}

void seissol::model::initializeSpecificLocalData( seissol::model::Material const& material,
                                                  seissol::model::LocalData* localData )
{
  auto sourceMatrix = init::ET::view::create(localData->sourceMatrix);
  getTransposedSourceCoefficientTensor(material, sourceMatrix);
}

void seissol::model::initializeSpecificNeighborData(  seissol::model::Material const&,
                                                      seissol::model::Material const (&)[4],
                                                      seissol::model::NeighborData* )
{
}
