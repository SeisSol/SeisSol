/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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
#ifndef ACOUSTIC_SETUP_H_
#define ACOUSTIC_SETUP_H_

#include "Model/Common.h"
#include "Kernels/Common.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"
#include "Numerical/Eigenvalues.h"

namespace seissol {
  namespace model {
    using Matrix44 = Eigen::Matrix<double, 4, 4>;

    template<typename T>
    inline void getTransposedCoefficientMatrix( AcousticMaterial const&  i_material,
                                                unsigned                i_dim,
                                                T&                      o_M )
      {
        o_M.setZero();

        real rhoInv = 1.0 / i_material.rho;

        switch (i_dim)
          {
            case 0:
              o_M(1,0) = i_material.lambda;
              o_M(0,1) = rhoInv;
              break;
        
            case 1:
              o_M(2,0) = i_material.lambda;
              o_M(0,2) = rhoInv;
              break;
        
            case 2:
              o_M(3,0) = i_material.lambda;
              o_M(0,3) = rhoInv;
              break;
              
            default:
              break;
          }
      }

    template<typename Tloc, typename Tneigh>
    inline void getTransposedGodunovState( AcousticMaterial const& local,
                                           AcousticMaterial const& neighbor,
                                           FaceType               faceType,
                                           Tloc&                  QgodLocal,
                                           Tneigh&                QgodNeighbor )
      {
         QgodNeighbor.setZero();

         // Eigenvectors are precomputed
         Matrix44 R = Matrix44::Zero();
         // scale for better condition number of R
         R(0,0) = std::sqrt(local.lambda * local.rho);
         R(1,0) = -local.lambda;
         R(0,1) = std::sqrt(neighbor.lambda * neighbor.rho);
         R(1,1) = neighbor.lambda;
         R(2,2) = local.lambda;
         R(3,3) = local.lambda;

         if (faceType == FaceType::FreeSurface) {
          //  MaterialType materialtype = testIfAcoustic(local.mu) ? MaterialType::Acoustic : MaterialType::Elastic;
          //  getTransposedFreeSurfaceGodunovState(materialtype, QgodLocal, QgodNeighbor, R);
           for (size_t i = 0; i < 4; i++) {
             for (size_t j = 0; j < 4; j++) {
               QgodNeighbor(i,j) = std::numeric_limits<double>::signaling_NaN();
             }
           }
           QgodLocal.setZero();
           QgodLocal(0, 1) = -1 * R(1,0) * 1/R(0,0);
           QgodLocal(1, 1) = 1.0;
         } else {
           Matrix44 chi = Matrix44::Zero();
           chi(0,0) = 1.0;
       
           const auto godunov = ((R*chi)*R.inverse()).eval();
       
           // QgodLocal = I - QgodNeighbor
           for (unsigned i = 0; i < godunov.cols(); ++i) {
             for (unsigned j = 0; j < godunov.rows(); ++j) {
               QgodLocal(i,j) = -godunov(j,i);
               QgodNeighbor(i,j) = godunov(j,i);
             }
           }
           for (unsigned idx = 0; idx < 4; ++idx) {
             QgodLocal(idx,idx) += 1.0;
           }
         }
      }
  }
}
#endif
