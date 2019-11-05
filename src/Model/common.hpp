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
 
#ifndef MODEL_COMMON_HPP_
#define MODEL_COMMON_HPP_

#include <Eigen/Eigen>
#include <cmath>
#include <Initializer/typedefs.hpp>
#include <generated_code/init.h>
#include <Geometry/MeshDefinition.h>
#include <Numerical_aux/Transformation.h>
#include <iostream>

namespace seissol {
  namespace model {
    using Matrix99 = Eigen::Matrix<real, 9, 9>;

    template<typename T>
    void getTransposedElasticCoefficientMatrix( ElasticMaterial const&          i_material,
                                                unsigned                        i_dim,
                                                T&                              o_M );

    template<typename Tloc, typename Tneigh>
    void getTransposedElasticGodunovState( Material const&                      local,
                                           Material const&                      neighbor,
                                           enum ::faceType                      faceType,
                                           Tloc&                                QgodLocal,
                                           Tneigh&                              QgodNeighbor );
    
    template<typename T>
    void getTransposedBoundaryGodunovState( T& QgodLocal,
                                            T& QgodNeighbor,
                                            Matrix99& R);

  }
}

template<typename T>
void seissol::model::getTransposedElasticCoefficientMatrix( seissol::model::ElasticMaterial const&  i_material,
                                                            unsigned                                i_dim,
                                                            T&                                      o_M )
{
  o_M.setZero();

  real lambda2mu = i_material.lambda + 2.0 * i_material.mu;
  real rhoInv = 1.0 / i_material.rho;

  switch (i_dim)
  {
    case 0:
      o_M(6,0) = -lambda2mu;
      o_M(6,1) = -i_material.lambda;
      o_M(6,2) = -i_material.lambda;
      o_M(7,3) = -i_material.mu;
      o_M(8,5) = -i_material.mu;
      o_M(0,6) = -rhoInv;
      o_M(3,7) = -rhoInv;
      o_M(5,8) = -rhoInv;
      break;

    case 1:
      o_M(7,0) = -i_material.lambda;
      o_M(7,1) = -lambda2mu;
      o_M(7,2) = -i_material.lambda;
      o_M(6,3) = -i_material.mu;
      o_M(8,4) = -i_material.mu;
      o_M(3,6) = -rhoInv;
      o_M(1,7) = -rhoInv;
      o_M(4,8) = -rhoInv;
      break;

    case 2:
      o_M(8,0) = -i_material.lambda;
      o_M(8,1) = -i_material.lambda;
      o_M(8,2) = -lambda2mu;
      o_M(7,4) = -i_material.mu;
      o_M(6,5) = -i_material.mu;
      o_M(5,6) = -rhoInv;
      o_M(4,7) = -rhoInv;
      o_M(2,8) = -rhoInv;
      break;
      
    default:
      break;
  }
}

template<typename T>
void seissol::model::getTransposedBoundaryGodunovState( T&                         QgodLocal,
                                                        T&                         QgodNeighbor,
                                                        Eigen::Matrix<real, 9, 9>& R)
{
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      QgodNeighbor(i,j) = std::numeric_limits<double>::signaling_NaN();
    }
  }

  QgodLocal.setZero();
  std::array<std::array<int, 2>, 3> traction_indices = {{{0,0}, {1,3}, {2,5}}};
  std::array<std::array<int, 2>, 3> velocity_indices = {{{0,6}, {1,7}, {2,8}}};
  //Eigen does not have fancy slicing options :(
  using Matrix33 = Eigen::Matrix<real, 3, 3>;
  Matrix33 R11 = Matrix33::Zero(); 
  for (auto &t: traction_indices) {
    for (int i = 0; i < 3; i++) {
      R11(t[0], i) = R(t[1],i);
    }
  }
  Matrix33 R21 = Matrix33::Zero();
  for (auto &v: velocity_indices) {
    for (int i = 0; i < 3; i++) {
      R21(v[0], i) = R(v[1],i);
    }
  }
  auto S = - (R21 * R11.inverse()).eval();

  //set lower left block
  for (auto &t: traction_indices) {
    for (auto &v: velocity_indices) {
      QgodLocal(t[1], v[1]) = S(v[0],t[0]);
    }
  }
  //set lower right block
  for (auto &v : velocity_indices) {
    QgodLocal(v[1], v[1]) = 1.0;
  }
}

template<typename Tloc, typename Tneigh>
void seissol::model::getTransposedElasticGodunovState( Material const&                      local,
    Material const&                      neighbor,
    enum ::faceType                      faceType,
    Tloc&                                QgodLocal,
    Tneigh&                              QgodNeighbor )
{
  QgodNeighbor.setZero();

  Matrix99 R = Matrix99::Zero();

  //eigenvectors have been precalculated
  R(0,0) = local.lambda + 2*local.mu;
  R(0,8) = neighbor.lambda + 2*neighbor.mu;
  R(1,0) = local.lambda;
  R(1,4) = 1;
  R(1,8) = neighbor.lambda;
  R(2,0) = local.lambda;
  R(2,5) = 1;
  R(2,8) = neighbor.lambda;
  R(3,1) = local.mu;
  R(3,7) = neighbor.mu;
  R(4,3) = 1;
  R(5,2) = local.mu;
  R(5,6) = neighbor.mu;
  R(6,0) = sqrt((local.lambda + 2*local.mu)/local.rho);
  R(6,8) = -sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho);
  R(7,1) = sqrt(local.mu/local.rho);
  R(7,7) = -sqrt(neighbor.mu/neighbor.rho);
  R(8,2) = sqrt(local.mu/local.rho);
  R(8,6) = -sqrt(neighbor.mu/neighbor.rho);

  if(faceType == freeSurface) {
    getTransposedBoundaryGodunovState(QgodLocal, QgodNeighbor, R);

  } else {
    Matrix99 R_inv = Matrix99::Zero();
    R_inv(0,0) = 1/(local.lambda + 2*local.mu) - (neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/((local.lambda + 2*local.mu)*(local.lambda + 2*local.mu)*((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho)));
    R_inv(0,6) = (neighbor.lambda + 2*neighbor.mu)/((local.lambda + 2*local.mu)*((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho)));
    R_inv(1,3) = 1/local.mu - neighbor.mu*sqrt(local.mu/local.rho)/(local.mu*local.mu*(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho)));
    R_inv(1,7) = neighbor.mu/(local.mu*(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho)));
    R_inv(2,5) = 1/local.mu - neighbor.mu*sqrt(local.mu/local.rho)/(local.mu*local.mu*(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho)));
    R_inv(2,8) = neighbor.mu/(local.mu*(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho)));
    R_inv(3,4) = 1;
    R_inv(4,0) = -local.lambda/(local.lambda + 2*local.mu) + (local.lambda*(neighbor.lambda + 2*neighbor.mu)/(local.lambda + 2*local.mu) - neighbor.lambda)*sqrt((local.lambda + 2*local.mu)/local.rho)/((local.lambda + 2*local.mu)*((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho)));
    R_inv(4,1) = 1;
    R_inv(4,6) = -(local.lambda*(neighbor.lambda + 2*neighbor.mu)/(local.lambda + 2*local.mu) - neighbor.lambda)/((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho));
    R_inv(5,0) = -local.lambda/(local.lambda + 2*local.mu) + (local.lambda*(neighbor.lambda + 2*neighbor.mu)/(local.lambda + 2*local.mu) - neighbor.lambda)*sqrt((local.lambda + 2*local.mu)/local.rho)/((local.lambda + 2*local.mu)*((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho)));
    R_inv(5,2) = 1;
    R_inv(5,6) = -(local.lambda*(neighbor.lambda + 2*neighbor.mu)/(local.lambda + 2*local.mu) - neighbor.lambda)/((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho));
    R_inv(6,5) = sqrt(local.mu/local.rho)/(local.mu*(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho)));
    R_inv(6,8) = -1/(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho));
    R_inv(7,3) = sqrt(local.mu/local.rho)/(local.mu*(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho)));
    R_inv(7,7) = -1/(neighbor.mu*sqrt(local.mu/local.rho)/local.mu + sqrt(neighbor.mu/neighbor.rho));
    R_inv(8,0) = sqrt((local.lambda + 2*local.mu)/local.rho)/((local.lambda + 2*local.mu)*((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho)));
    R_inv(8,6) = -1/((neighbor.lambda + 2*neighbor.mu)*sqrt((local.lambda + 2*local.mu)/local.rho)/(local.lambda + 2*local.mu) + sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho));

    Matrix99 chi = Matrix99::Zero();
    chi(0,0) = 1.0;
    chi(1,1) = 1.0;
    chi(2,2) = 1.0;

    auto godunov = ((R*chi)*R_inv).eval();

    // QgodLocal = I - QgodNeighbor
    for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
      for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
        QgodLocal(i,j) = -godunov(j,i);
        QgodNeighbor(i,j) = godunov(j,i);
      }
    }  
    for (unsigned idx = 0; idx < 9; ++idx) {
      QgodLocal(idx,idx) += 1.0;
    }
  }
}

#endif
