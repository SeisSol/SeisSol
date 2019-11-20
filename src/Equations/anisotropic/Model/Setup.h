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

#ifndef ANISOTROPIC_SETUP_H_
#define ANISOTROPIC_SETUP_H_

#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>

#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>
#include <generated_code/init.h>

namespace seissol {
  namespace model {

    template<typename T>
      inline void getTransposedCoefficientMatrix( AnisotropicMaterial const&  i_material,
          unsigned                    i_dim,
          T&                          o_M )
      {
        o_M.setZero();

        real rhoInv = 1.0 / i_material.rho;

        switch (i_dim)
        {
          case 0:
            o_M(6,0) = -i_material.c11;
            o_M(7,0) = -i_material.c16;
            o_M(8,0) = -i_material.c15;
            o_M(6,1) = -i_material.c12;
            o_M(7,1) = -i_material.c26;
            o_M(8,1) = -i_material.c25;
            o_M(6,2) = -i_material.c13;
            o_M(7,2) = -i_material.c36;
            o_M(8,2) = -i_material.c35;
            o_M(6,3) = -i_material.c16;
            o_M(7,3) = -i_material.c66;
            o_M(8,3) = -i_material.c56;
            o_M(6,4) = -i_material.c14;
            o_M(7,4) = -i_material.c46;
            o_M(8,4) = -i_material.c45;
            o_M(6,5) = -i_material.c15;
            o_M(7,5) = -i_material.c56;
            o_M(8,5) = -i_material.c55;
            o_M(0,6) = -rhoInv;
            o_M(3,7) = -rhoInv;
            o_M(5,8) = -rhoInv;
            break;

          case 1:
            o_M(6,0) = -i_material.c16;
            o_M(7,0) = -i_material.c12;
            o_M(8,0) = -i_material.c14;
            o_M(6,1) = -i_material.c26;
            o_M(7,1) = -i_material.c22;
            o_M(8,1) = -i_material.c24;
            o_M(6,2) = -i_material.c36;
            o_M(7,2) = -i_material.c23;
            o_M(8,2) = -i_material.c34;
            o_M(6,3) = -i_material.c66;
            o_M(7,3) = -i_material.c26;
            o_M(8,3) = -i_material.c46;
            o_M(6,4) = -i_material.c46;
            o_M(7,4) = -i_material.c24;
            o_M(8,4) = -i_material.c44;
            o_M(6,5) = -i_material.c56;
            o_M(7,5) = -i_material.c25;
            o_M(8,5) = -i_material.c45;
            o_M(3,6) = -rhoInv;
            o_M(1,7) = -rhoInv;
            o_M(4,8) = -rhoInv;
            break;

          case 2:
            o_M(6,0) = -i_material.c15;
            o_M(7,0) = -i_material.c14;
            o_M(8,0) = -i_material.c13;
            o_M(6,1) = -i_material.c25;
            o_M(7,1) = -i_material.c24;
            o_M(8,1) = -i_material.c23;
            o_M(6,2) = -i_material.c35;
            o_M(7,2) = -i_material.c34;
            o_M(8,2) = -i_material.c33;
            o_M(6,3) = -i_material.c56;
            o_M(7,3) = -i_material.c46;
            o_M(8,3) = -i_material.c36;
            o_M(6,4) = -i_material.c45;
            o_M(7,4) = -i_material.c44;
            o_M(8,4) = -i_material.c34;
            o_M(6,5) = -i_material.c55;
            o_M(7,5) = -i_material.c45;
            o_M(8,5) = -i_material.c35;
            o_M(5,6) = -rhoInv;
            o_M(4,7) = -rhoInv;
            o_M(2,8) = -rhoInv;
            break;

          default:
            break;
        }
      }

    inline void getEigenBasisForAnisotropicMaterial( AnisotropicMaterial const&  local,
        AnisotropicMaterial const&  neighbor,
        Eigen::Matrix<real, 9, 9>&  R)
    {
      using Matrix33 = Eigen::Matrix<real, 3, 3, Eigen::ColMajor>;
      using Matrix63 = Eigen::Matrix<real, 6, 3, Eigen::ColMajor>;

      //Calculate Eigenvectors and Eigenvalues
      Eigen::SelfAdjointEigenSolver<Matrix33> saes;

      real aL[9];
      aL[0] = local.c11 / local.rho;  
      aL[1] = local.c16 / local.rho;  
      aL[2] = local.c15 / local.rho;  
      aL[3] = local.c16 / local.rho;  
      aL[4] = local.c66 / local.rho;  
      aL[5] = local.c56 / local.rho;  
      aL[6] = local.c15 / local.rho;  
      aL[7] = local.c56 / local.rho;  
      aL[8] = local.c55 / local.rho;  
      Matrix33 AL(aL);
      saes.compute(AL);
      auto eigenvaluesL = saes.eigenvalues();
      auto eigenvectorsL = saes.eigenvectors();

      real aN[9];
      aN[0] = neighbor.c11 / neighbor.rho;  
      aN[1] = neighbor.c16 / neighbor.rho;  
      aN[2] = neighbor.c15 / neighbor.rho;  
      aN[3] = neighbor.c16 / neighbor.rho;  
      aN[4] = neighbor.c66 / neighbor.rho;  
      aN[5] = neighbor.c56 / neighbor.rho;  
      aN[6] = neighbor.c15 / neighbor.rho;  
      aN[7] = neighbor.c56 / neighbor.rho;  
      aN[8] = neighbor.c55 / neighbor.rho;  
      Matrix33 AN(aN);
      saes.compute(AN);
      auto eigenvaluesN = saes.eigenvalues();
      auto eigenvectorsN = saes.eigenvectors();

      real a1L[18];
      a1L[0] = -local.c11;
      a1L[1] = -local.c12;
      a1L[2] = -local.c13;
      a1L[3] = -local.c16;
      a1L[4] = -local.c14;
      a1L[5] = -local.c15;
      a1L[6] = -local.c16;
      a1L[7] = -local.c26;
      a1L[8] = -local.c36;
      a1L[9] = -local.c66;
      a1L[10] = -local.c46;
      a1L[11] = -local.c56;
      a1L[12] = -local.c15;
      a1L[13] = -local.c25;
      a1L[14] = -local.c35;
      a1L[15] = -local.c46;
      a1L[16] = -local.c45;
      a1L[17] = -local.c55;
      Matrix63 A1L(a1L);

      real a1N[18];
      a1N[0] = -neighbor.c11;
      a1N[1] = -neighbor.c12;
      a1N[2] = -neighbor.c13;
      a1N[3] = -neighbor.c16;
      a1N[4] = -neighbor.c14;
      a1N[5] = -neighbor.c15;
      a1N[6] = -neighbor.c16;
      a1N[7] = -neighbor.c26;
      a1N[8] = -neighbor.c36;
      a1N[9] = -neighbor.c66;
      a1N[10] = -neighbor.c46;
      a1N[11] = -neighbor.c56;
      a1N[12] = -neighbor.c15;
      a1N[13] = -neighbor.c25;
      a1N[14] = -neighbor.c35;
      a1N[15] = -neighbor.c46;
      a1N[16] = -neighbor.c45;
      a1N[17] = -neighbor.c55;
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

      R <<  A1L * eigenvectorsL,   E1, A1N * eigenvectorsN, 
        -eigenvectorsL*lambdaL, E2, eigenvectorsN*lambdaN;
    }

    template<>
      inline void getTransposedGodunovState( AnisotropicMaterial const&       local,
          AnisotropicMaterial const&       neighbor,
          enum ::faceType                  faceType,
          init::QgodLocal::view::type&     QgodLocal,
          init::QgodNeighbor::view::type&  QgodNeighbor )
      {

        Matrix99 R = Matrix99::Zero();
        getEigenBasisForAnisotropicMaterial(local, neighbor, R);

        if(faceType == freeSurface) {
          getTransposedFreeSurfaceGodunovState(QgodLocal, QgodNeighbor, R);

        } else {
          Matrix99 chi = Matrix99::Zero();
          chi(0,0) = 1.0;
          chi(1,1) = 1.0;
          chi(2,2) = 1.0;

          const auto godunov = ((R*chi)*R.inverse()).eval();

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
        }
      }


    template<>
      inline void initializeSpecificLocalData( AnisotropicMaterial const&,
          LocalData* )
      {
      }

    template<>
      inline void initializeSpecificNeighborData(  AnisotropicMaterial const&,
          NeighborData* )
      {
      }

    template<>
      inline void setMaterial( double*               i_materialVal,
          int                   i_numMaterialVals,
          AnisotropicMaterial*  o_material )
      {
        assert(i_numMaterialVals == 22);
        o_material->rho = i_materialVal[0];
        o_material->c11 = i_materialVal[1];
        o_material->c12 = i_materialVal[2];
        o_material->c13 = i_materialVal[3];
        o_material->c14 = i_materialVal[4];
        o_material->c15 = i_materialVal[5];
        o_material->c16 = i_materialVal[6];
        o_material->c22 = i_materialVal[7];
        o_material->c23 = i_materialVal[8];
        o_material->c24 = i_materialVal[9];
        o_material->c25 = i_materialVal[10];
        o_material->c26 = i_materialVal[11];
        o_material->c33 = i_materialVal[12];
        o_material->c34 = i_materialVal[13];
        o_material->c35 = i_materialVal[14];
        o_material->c36 = i_materialVal[15];
        o_material->c44 = i_materialVal[16];
        o_material->c45 = i_materialVal[17];
        o_material->c46 = i_materialVal[18];
        o_material->c55 = i_materialVal[19];
        o_material->c56 = i_materialVal[20];
        o_material->c66 = i_materialVal[21];
      }

    inline AnisotropicMaterial getRotatedMaterialCoefficients(real i_N[36], AnisotropicMaterial& i_material) {
        AnisotropicMaterial o_material;
        o_material.rho = i_material.rho;
        using Matrix66 = Eigen::Matrix<real, 6, 6>;
        Matrix66 N = Matrix66(i_N);
        Matrix66 C = Matrix66();
        C(0,0) = i_material.c11;
        C(0,1) = i_material.c12;
        C(0,2) = i_material.c13;
        C(0,3) = i_material.c14;
        C(0,4) = i_material.c15;
        C(0,5) = i_material.c16;
        C(1,0) = i_material.c12;
        C(1,1) = i_material.c22;
        C(1,2) = i_material.c23;
        C(1,3) = i_material.c24;
        C(1,4) = i_material.c25;
        C(1,5) = i_material.c26;
        C(2,0) = i_material.c13;
        C(2,1) = i_material.c23;
        C(2,2) = i_material.c33;
        C(2,3) = i_material.c34;
        C(2,4) = i_material.c35;
        C(2,5) = i_material.c36;
        C(3,0) = i_material.c14;
        C(3,1) = i_material.c24;
        C(3,2) = i_material.c34;
        C(3,3) = i_material.c44;
        C(3,4) = i_material.c45;
        C(3,5) = i_material.c46;
        C(4,0) = i_material.c15;
        C(4,1) = i_material.c25;
        C(4,2) = i_material.c35;
        C(4,3) = i_material.c45;
        C(4,4) = i_material.c55;
        C(4,5) = i_material.c56;
        C(5,0) = i_material.c16;
        C(5,1) = i_material.c26;
        C(5,2) = i_material.c36;
        C(5,3) = i_material.c46;
        C(5,4) = i_material.c56;
        C(5,5) = i_material.c66;

        Matrix66 rotatedC = N.transpose()*C*N;

        o_material.c11 = rotatedC(0,0);
        o_material.c12 = rotatedC(0,1);
        o_material.c13 = rotatedC(0,2);
        o_material.c14 = rotatedC(0,3);
        o_material.c15 = rotatedC(0,4);
        o_material.c16 = rotatedC(0,5);
        o_material.c22 = rotatedC(1,1);
        o_material.c23 = rotatedC(1,2);
        o_material.c24 = rotatedC(1,3);
        o_material.c25 = rotatedC(1,4);
        o_material.c26 = rotatedC(1,5);
        o_material.c33 = rotatedC(2,2);
        o_material.c34 = rotatedC(2,3);
        o_material.c35 = rotatedC(2,4);
        o_material.c36 = rotatedC(2,5);
        o_material.c44 = rotatedC(3,3);
        o_material.c45 = rotatedC(3,4);
        o_material.c46 = rotatedC(3,5);
        o_material.c55 = rotatedC(4,4);
        o_material.c56 = rotatedC(4,5);
        o_material.c66 = rotatedC(5,5);
        return o_material;
      }

  }
}

#endif
