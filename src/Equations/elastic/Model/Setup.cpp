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

#include <generated_code/init.h>

void seissol::model::getTransposedCoefficientMatrix( Material const&                i_material,
                                                     unsigned                       i_dim,
                                                     init::star::view<0>::type&     AT )
{
  seissol::model::getTransposedElasticCoefficientMatrix(i_material, i_dim, AT);
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
  assert(i_numMaterialVals == 3);
 
  o_material->rho = i_materialVal[0];
  o_material->mu = i_materialVal[1];
  o_material->lambda = i_materialVal[2];
}


void seissol::model::getFaceRotationMatrix( VrtxCoords const i_normal,
                                            VrtxCoords const i_tangent1,
                                            VrtxCoords const i_tangent2,
                                            init::T::view::type&    o_T,
                                            init::Tinv::view::type& o_Tinv )
{
  o_T.setZero();
  o_Tinv.setZero();
  
  unsigned origin0[]{0, 0};
  unsigned shape6[]{6, 6};
  auto TBlock0 = o_T.block(origin0, shape6);
  auto TinvBlock0 = o_Tinv.block(origin0, shape6);
  
  unsigned origin6[]{0, 0};
  unsigned shape3[]{3, 3};
  auto TBlock1 = o_T.block(origin6, shape3);
  auto TinvBlock1 = o_Tinv.block(origin6, shape3);
  
  seissol::transformations::symmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, TBlock0);
  seissol::transformations::tensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, TBlock1);
  
  seissol::transformations::inverseSymmetricTensor2RotationMatrix(i_normal, i_tangent1, i_tangent2, TinvBlock0);
  seissol::transformations::inverseTensor1RotationMatrix(i_normal, i_tangent1, i_tangent2, TinvBlock1);
}

void seissol::model::initializeSpecificLocalData( seissol::model::Material const&,
                                                  seissol::model::LocalData* )
{
}

void seissol::model::initializeSpecificNeighborData(  seissol::model::Material const&,
                                                      seissol::model::Material const (&)[4],
                                                      seissol::model::NeighborData* )
{
}
