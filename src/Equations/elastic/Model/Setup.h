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
#ifndef ELASTIC_SETUP_H_
#define ELASTIC_SETUP_H_

#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>
#include <generated_code/init.h>

namespace seissol {
  namespace model {

    template<typename T>
      inline void getTransposedCoefficientMatrix( ElasticMaterial const&  i_material,
                                                  unsigned                i_dim,
                                                  T&                      o_M )
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

    template<typename Tloc, typename Tneigh>
      inline void getTransposedGodunovState( ElasticMaterial const& local,
                                             ElasticMaterial const& neighbor,
                                             enum ::faceType        faceType,
                                             Tloc&                  QgodLocal,
                                             Tneigh&                QgodNeighbor )
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
          getTransposedFreeSurfaceGodunovState(QgodLocal, QgodNeighbor, R);

        } else {
          Matrix99 R_inv = Matrix99::Zero();
          //we do not need to compute all of R_inv as we will multiply it with the indicator chi later
          //which extracts the first three rows of R_inv

          //We can exploit that R only couples 2 values to each other to easily compute an analytic solution
          //inv_xy computes the (x,y)th entry of ((local_kappa, neighbor_kappa),(local_c, -neighbor_c))^-1
          auto inv_00 = [](real local_kappa, real neighbor_kappa, real local_c, real neighbor_c) {
            return neighbor_c / (local_kappa * neighbor_c + neighbor_kappa * local_c);  
          };
          auto inv_01 = [](real local_kappa, real neighbor_kappa, real local_c, real neighbor_c) {
            return neighbor_kappa / (local_kappa * neighbor_c + neighbor_kappa * local_c);  
          };


          R_inv(0,0) = inv_00(local.lambda + 2*local.mu, neighbor.lambda + 2*neighbor.mu, sqrt((local.lambda + 2*local.mu)/local.rho), sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho));
          R_inv(0,6) = inv_01(local.lambda + 2*local.mu, neighbor.lambda + 2*neighbor.mu, sqrt((local.lambda + 2*local.mu)/local.rho), sqrt((neighbor.lambda + 2*neighbor.mu)/neighbor.rho));
          R_inv(1,3) = inv_00(local.mu, neighbor.mu, sqrt(local.mu/local.rho), sqrt(neighbor.mu/neighbor.rho));
          R_inv(1,7) = inv_01(local.mu, neighbor.mu, sqrt(local.mu/local.rho), sqrt(neighbor.mu/neighbor.rho));
          R_inv(2,5) = inv_00(local.mu, neighbor.mu, sqrt(local.mu/local.rho), sqrt(neighbor.mu/neighbor.rho));
          R_inv(2,8) = inv_01(local.mu, neighbor.mu, sqrt(local.mu/local.rho), sqrt(neighbor.mu/neighbor.rho));

          Matrix99 chi = Matrix99::Zero();
          chi(0,0) = 1.0;
          chi(1,1) = 1.0;
          chi(2,2) = 1.0;

          const auto godunov = ((R*chi)*R_inv).eval();

          // QgodLocal = I - QgodNeighbor
          for (unsigned i = 0; i < godunov.cols(); ++i) {
            for (unsigned j = 0; j < godunov.rows(); ++j) {
              QgodLocal(i,j) = -godunov(j,i);
              QgodNeighbor(i,j) = godunov(j,i);
            }
          }  
          for (unsigned idx = 0; idx < 9; ++idx) {
            QgodLocal(idx,idx) += 1.0;
          }
        }
      }

    template<>
    inline void initializeSpecificLocalData( ElasticMaterial const&,
                                      LocalData* )
    {
    }

    template<>
    inline void initializeSpecificNeighborData( ElasticMaterial const&,
                                         NeighborData* )
    {
    }

    template<>
    inline void setMaterial( double*           i_materialVal,
                      int               i_numMaterialVals,
                      ElasticMaterial*  o_material )
    {
      assert(i_numMaterialVals == 3);

      o_material->rho = i_materialVal[0];
      o_material->mu = i_materialVal[1];
      o_material->lambda = i_materialVal[2];
    }

  }
}
#endif
