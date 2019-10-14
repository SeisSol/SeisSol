/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
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
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>

#include <Model/Setup.h>
#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>
#include <generated_code/init.h>
#include <iostream>
template<typename T>
void getTransposedAnisotropicCoefficientMatrix( seissol::model::Material const&  i_material,
                                                unsigned                         i_dim,
                                                T&                               o_M )
{
  o_M.setZero();
  
  real rhoInv = 1.0 / i_material.rho;

  switch (i_dim)
  {
    case 0:
      o_M(6,0) = -i_material.c(0,0);
      o_M(7,0) = -i_material.c(0,5);
      o_M(8,0) = -i_material.c(0,4);
      o_M(6,1) = -i_material.c(0,1);
      o_M(7,1) = -i_material.c(1,5);
      o_M(8,1) = -i_material.c(1,4);
      o_M(6,2) = -i_material.c(0,2);
      o_M(7,2) = -i_material.c(2,5);
      o_M(8,2) = -i_material.c(2,4);
      o_M(6,3) = -i_material.c(0,5);
      o_M(7,3) = -i_material.c(5,5);
      o_M(8,3) = -i_material.c(4,5);
      o_M(6,4) = -i_material.c(0,3);
      o_M(7,4) = -i_material.c(3,5);
      o_M(8,4) = -i_material.c(3,4);
      o_M(6,5) = -i_material.c(0,4);
      o_M(7,5) = -i_material.c(4,5);
      o_M(8,5) = -i_material.c(4,4);
      o_M(0,6) = -rhoInv;
      o_M(3,7) = -rhoInv;
      o_M(5,8) = -rhoInv;
      break;

    case 1:
      o_M(6,0) = -i_material.c(0,5);
      o_M(7,0) = -i_material.c(0,1);
      o_M(8,0) = -i_material.c(0,3);
      o_M(6,1) = -i_material.c(1,5);
      o_M(7,1) = -i_material.c(1,1);
      o_M(8,1) = -i_material.c(1,3);
      o_M(6,2) = -i_material.c(2,5);
      o_M(7,2) = -i_material.c(1,2);
      o_M(8,2) = -i_material.c(2,3);
      o_M(6,3) = -i_material.c(5,5);
      o_M(7,3) = -i_material.c(1,5);
      o_M(8,3) = -i_material.c(3,5);
      o_M(6,4) = -i_material.c(3,5);
      o_M(7,4) = -i_material.c(1,3);
      o_M(8,4) = -i_material.c(3,3);
      o_M(6,5) = -i_material.c(4,5);
      o_M(7,5) = -i_material.c(1,4);
      o_M(8,5) = -i_material.c(3,4);
      o_M(3,6) = -rhoInv;
      o_M(1,7) = -rhoInv;
      o_M(4,8) = -rhoInv;
      break;

    case 2:
      o_M(6,0) = -i_material.c(0,4);
      o_M(7,0) = -i_material.c(0,3);
      o_M(8,0) = -i_material.c(0,2);
      o_M(6,1) = -i_material.c(1,4);
      o_M(7,1) = -i_material.c(1,3);
      o_M(8,1) = -i_material.c(1,2);
      o_M(6,2) = -i_material.c(2,4);
      o_M(7,2) = -i_material.c(2,3);
      o_M(8,2) = -i_material.c(2,2);
      o_M(6,3) = -i_material.c(4,5);
      o_M(7,3) = -i_material.c(3,5);
      o_M(8,3) = -i_material.c(2,5);
      o_M(6,4) = -i_material.c(3,4);
      o_M(7,4) = -i_material.c(3,3);
      o_M(8,4) = -i_material.c(2,3);
      o_M(6,5) = -i_material.c(4,4);
      o_M(7,5) = -i_material.c(3,4);
      o_M(8,5) = -i_material.c(2,4);
      o_M(5,6) = -rhoInv;
      o_M(4,7) = -rhoInv;
      o_M(2,8) = -rhoInv;
      break;
      
    default:
      break;
  }
}

void seissol::model::getTransposedGodunovState( Material const&                   local,
                                                Material const&                   neighbor,
                                                enum ::faceType                   faceType,
                                                init::QgodLocal::view::type&      QgodLocal,
                                                init::QgodNeighbor::view::type&   QgodNeighbor )
{
  QgodNeighbor.setZero();
  
  using Matrix33 = Eigen::Matrix<real, 3, 3, Eigen::ColMajor>;
  using Matrix63 = Eigen::Matrix<real, 6, 3, Eigen::ColMajor>;
  using Matrix99 = Eigen::Matrix<real, 9, 9, Eigen::ColMajor>;
  
//Calculate Eigenvectors and Eigenvalues
  Eigen::SelfAdjointEigenSolver<Matrix33> saes;
  
  real aL[9];
  aL[0] = local.c(0,0) / local.rho;  
  aL[1] = local.c(0,5) / local.rho;  
  aL[2] = local.c(0,4) / local.rho;  
  aL[3] = local.c(0,5) / local.rho;  
  aL[4] = local.c(5,5) / local.rho;  
  aL[5] = local.c(4,5) / local.rho;  
  aL[6] = local.c(0,4) / local.rho;  
  aL[7] = local.c(4,5) / local.rho;  
  aL[8] = local.c(4,4) / local.rho;  
  Matrix33 AL(aL);
  saes.compute(AL);
  auto eigenvaluesL = saes.eigenvalues();
  auto eigenvectorsL = saes.eigenvectors();

  real aN[9];
  aN[0] = neighbor.c(0,0) / neighbor.rho;  
  aN[1] = neighbor.c(0,5) / neighbor.rho;  
  aN[2] = neighbor.c(0,4) / neighbor.rho;  
  aN[3] = neighbor.c(0,5) / neighbor.rho;  
  aN[4] = neighbor.c(5,5) / neighbor.rho;  
  aN[5] = neighbor.c(4,5) / neighbor.rho;  
  aN[6] = neighbor.c(0,4) / neighbor.rho;  
  aN[7] = neighbor.c(4,5) / neighbor.rho;  
  aN[8] = neighbor.c(4,4) / neighbor.rho;  
  Matrix33 AN(aN);
  saes.compute(AN);
  auto eigenvaluesN = saes.eigenvalues();
  auto eigenvectorsN = saes.eigenvectors();

  real a1L[18];
  a1L[0] = -local.c(0,0);
  a1L[1] = -local.c(0,1);
  a1L[2] = -local.c(0,2);
  a1L[3] = -local.c(0,5);
  a1L[4] = -local.c(0,3);
  a1L[5] = -local.c(0,4);
  a1L[6] = -local.c(0,5);
  a1L[7] = -local.c(1,5);
  a1L[8] = -local.c(2,5);
  a1L[9] = -local.c(5,5);
  a1L[10] = -local.c(3,5);
  a1L[11] = -local.c(4,5);
  a1L[12] = -local.c(0,4);
  a1L[13] = -local.c(1,4);
  a1L[14] = -local.c(2,4);
  a1L[15] = -local.c(3,5);
  a1L[16] = -local.c(3,4);
  a1L[17] = -local.c(4,4);
  Matrix63 A1L(a1L);
  
  real a1N[18];
  a1N[0] = -neighbor.c(0,0);
  a1N[1] = -neighbor.c(0,1);
  a1N[2] = -neighbor.c(0,2);
  a1N[3] = -neighbor.c(0,5);
  a1N[4] = -neighbor.c(0,3);
  a1N[5] = -neighbor.c(0,4);
  a1N[6] = -neighbor.c(0,5);
  a1N[7] = -neighbor.c(1,5);
  a1N[8] = -neighbor.c(2,5);
  a1N[9] = -neighbor.c(5,5);
  a1N[10] = -neighbor.c(3,5);
  a1N[11] = -neighbor.c(4,5);
  a1N[12] = -neighbor.c(0,4);
  a1N[13] = -neighbor.c(1,4);
  a1N[14] = -neighbor.c(2,4);
  a1N[15] = -neighbor.c(3,5);
  a1N[16] = -neighbor.c(3,4);
  a1N[17] = -neighbor.c(4,4);
  Matrix63 A1N(a1N);
  
  for(unsigned i = 0; i < 3; i++) {
    eigenvaluesL(i) = sqrt(eigenvaluesL(i));
  }
  auto lambdaL = Eigen::Matrix<real, 3, 3, Eigen::ColMajor>(eigenvaluesL.asDiagonal());
  for(unsigned i = 0; i < 3; i++) {
     eigenvaluesN(i) = sqrt(eigenvaluesN(i));
  }
  auto lambdaN = Matrix33(eigenvaluesN.asDiagonal());

  Matrix63 E1 = Matrix63::Zero();
  E1(1,0) = 1;
  E1(2,1) = 1;
  E1(4,2) = 1;
  Matrix33 E2 = Matrix33::Zero();
  
  Matrix99 R;
  R <<  A1L * eigenvectorsL,   E1, A1N * eigenvectorsN, 
       -eigenvectorsL*lambdaL, E2, eigenvectorsN*lambdaN;
  Eigen::Matrix<real, 9, 1> diag = Eigen::Matrix<real, 9, 1>::Zero();
  diag(0) = 1.0;
  diag(1) = 1.0;
  diag(2) = 1.0;
  auto chi = Matrix99(diag.asDiagonal());
  Matrix99 godunov = ((R*chi)*R.inverse()).eval();
  
  // QgodLocal = I - QgodNeighbor
  for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
    for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
      QgodLocal(i,j) = -godunov(j,i);
      QgodNeighbor(i,j) = godunov(j,i);
    }
  }  
  for (unsigned idx = 0; idx < QgodLocal.shape(0) && idx < QgodLocal.shape(1); ++idx) {
    QgodLocal(idx,idx) += 1.0;
  }
  applyBoundaryConditionToElasticFluxSolver(faceType, QgodNeighbor);
}

void seissol::model::getPlaneWaveOperator(  Material const& material,
                                            double const n[3],
                                            std::complex<real> Mdata[9 * 9] )
{
  yateto::DenseTensorView<2,std::complex<real>> M(Mdata, {9, 9});
  M.setZero();

  real data[9 * 9];
  yateto::DenseTensorView<2,real> Coeff(data, {9, 9});

  for (unsigned d = 0; d < 3; ++d) {
    Coeff.setZero();
    getTransposedAnisotropicCoefficientMatrix(material, d, Coeff);

    for (unsigned i = 0; i < 9; ++i) {
      for (unsigned j = 0; j < 9; ++j) {
        M(i,j) += n[d] * Coeff(j,i);
      }
    }
  }
}

void seissol::model::getTransposedCoefficientMatrix( Material const&                i_material,
                                                     unsigned                       i_dim,
                                                     init::star::view<0>::type&     AT )
{
  getTransposedAnisotropicCoefficientMatrix( i_material, i_dim, AT);
}
 
void seissol::model::setMaterial( double* i_materialVal,
                                  int i_numMaterialVals,
                                  seissol::model::Material* o_material )
{
  assert(i_numMaterialVals == 22);
  o_material->rho = i_materialVal[0];
  std::copy(i_materialVal+1, i_materialVal+21, o_material->c_store);
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

void seissol::model::initializeSpecificLocalData( seissol::model::Material const&,
                                                  seissol::model::LocalData* )
{
}

void seissol::model::initializeSpecificNeighborData(  seissol::model::Material const&,
                                                      seissol::model::Material const (&)[4],
                                                      seissol::model::NeighborData* )
{
}
