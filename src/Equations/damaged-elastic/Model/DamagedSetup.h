#ifndef DAMAGED_SETUP_H_
#define DAMAGED_SETUP_H_

#include <Model/common.hpp>
#include <Kernels/common.hpp>
#include <Numerical_aux/Transformation.h>
#include <generated_code/init.h>

#include <cmath>

namespace seissol {
  namespace model {
    using Matrix1010 = Eigen::Matrix<double, 11, 11>;

    template<typename T>
    inline void getTransposedCoefficientMatrix( DamagedElasticMaterial  const&  i_material,
                                                unsigned                i_dim,
                                                T&                      o_M )
      {
        o_M.setZero();

        // real epsInit = -1e-1; // eps_xx0
        real lambda2muInvRho = (i_material.lambda + 2.0 * i_material.mu)/i_material.rho;
        real lambdaInvRho = i_material.lambda / i_material.rho;
        real muInvRho = i_material.mu / i_material.rho;

        real I1 = i_material.epsxx_alpha + i_material.epsyy_alpha + i_material.epszz_alpha;
        real I2 = i_material.epsxx_alpha*i_material.epsxx_alpha
          + i_material.epsyy_alpha*i_material.epsyy_alpha
          + i_material.epszz_alpha*i_material.epszz_alpha
          + 2*i_material.epsxy_alpha*i_material.epsxy_alpha
          + 2*i_material.epsyz_alpha*i_material.epsyz_alpha
          + 2*i_material.epszx_alpha*i_material.epszx_alpha;

        // real xi = I1 / std::sqrt(I2 + 1e-20);
        // real xiInv = 1 / (xi+1e-2);

        real xi;
          if (I2 > 1e-30){
            xi = I1 / std::sqrt(I2);
          } else{
            xi = 0.0;
          }

          real xiInv;
          if ( std::abs(xi) > 1e-1){
            xiInv = 1 / xi;
          } else{
            xiInv = 0.0;
          }

        /*
        For strain-vel + damage formula
        {eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_zz, vx, vy, vz, alpha}
        */
        switch (i_dim)
          {
            case 0:
              o_M(6,0) = -1.0;
              o_M(7,3) = -0.5;
              o_M(8,5) = -0.5;
              o_M(0,6) = -lambda2muInvRho;
              o_M(1,6) = -lambdaInvRho;
              o_M(2,6) = -lambdaInvRho;
              if (!testIfAcoustic(i_material.mu)) {
                o_M(3,7) = -2.0*muInvRho;
                o_M(5,8) = -2.0*muInvRho;
              }
              o_M(9,6) = (i_material.gammaR*std::sqrt(I2)
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epsxx_alpha )
                  / i_material.rho;
              o_M(9,7) = (0
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epsxy_alpha )
                  / i_material.rho;
              o_M(9,8) = (0
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epszx_alpha )
                  / i_material.rho;
              break;

            case 1:
              o_M(7,1) = -1.0;
              o_M(6,3) = -0.5;
              o_M(8,4) = -0.5;
              o_M(1,7) = -lambda2muInvRho;
              o_M(0,7) = -lambdaInvRho;
              o_M(2,7) = -lambdaInvRho;
              if (!testIfAcoustic(i_material.mu)) {
                o_M(3,6) = -2.0*muInvRho;
                o_M(4,8) = -2.0*muInvRho;
              }
              o_M(9,6) = (0
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epsxy_alpha )
                  / i_material.rho;
              o_M(9,7) = (i_material.gammaR*std::sqrt(I2)
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epsyy_alpha )
                  / i_material.rho;
              o_M(9,8) = (0
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epsyz_alpha )
                  / i_material.rho;
              break;

            case 2:
              o_M(8,2) = -1.0;
              o_M(7,4) = -0.5;
              o_M(6,5) = -0.5;
              o_M(2,8) = -lambda2muInvRho;
              o_M(0,8) = -lambdaInvRho;
              o_M(1,8) = -lambdaInvRho;
              if (!testIfAcoustic(i_material.mu)) {
                o_M(5,6) = -2.0*muInvRho;
                o_M(4,7) = -2.0*muInvRho;
              }
              o_M(9,6) = (0
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epszx_alpha )
                  / i_material.rho;
              o_M(9,7) = (0
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epsyz_alpha )
                  / i_material.rho;
              o_M(9,8) = (i_material.gammaR*std::sqrt(I2)
                + i_material.gammaR*(xi+2.0*i_material.xi0)*i_material.epszz_alpha )
                  / i_material.rho;
              break;

            default:
              break;
          }
      }

    template<typename Tloc, typename Tneigh>
    inline void getTransposedGodunovState( DamagedElasticMaterial const& local,
                                           DamagedElasticMaterial const& neighbor,
                                           FaceType               faceType,
                                           Tloc&                  QgodLocal,
                                           Tneigh&                QgodNeighbor )
      {
         QgodNeighbor.setZero();

         // Eigenvectors are precomputed
         Matrix1010 R = Matrix1010::Zero();

        /*
        For strain-vel + damage formula
        {-cp, -cs, -cs, 0, 0, 0, cs, cs, cp, 0}
        */
        if (testIfAcoustic(local.mu)) {
           R(0,0) = local.lambda;
           R(1,0) = local.lambda;
           R(2,0) = local.lambda;
           R(6,0) = std::sqrt((local.lambda) / local.rho);

           // scale for better condition number of R
           R(3,1) = local.lambda;
           R(5,2) = local.lambda;
         } else {
           R(0,0) = 1;
           R(6,0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);

           R(5,1) = 0.5;
           R(8,1) = std::sqrt(local.mu / local.rho);

           R(3,2) = 0.5;
           R(7,2) = std::sqrt(local.mu / local.rho);
         }

         // scale for better condition number of R
         R(4,3) = 1;

         R(0,4) = -1;
         R(2,4) = (local.lambda + 2*local.mu)/local.lambda;

         R(0,5) = -1;
         R(1,5) = (local.lambda + 2*local.mu)/local.lambda;

         if (testIfAcoustic(neighbor.mu)) {
           // scale for better condition number of R
           R(7,6) = neighbor.lambda;
           R(8,7) = neighbor.lambda;

           R(0,8) = neighbor.lambda;
           R(1,8) = neighbor.lambda;
           R(2,8) = neighbor.lambda;
           R(6,8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
         } else {
           R(5,6) = 0.5;
           R(8,6) = -std::sqrt(neighbor.mu / neighbor.rho);

           R(3,7) = 0.5;
           R(7,7) = -std::sqrt(neighbor.mu / neighbor.rho);

           R(0,8) = 1.0;
           R(6,8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);
         }
         #if USE_DAMAGEDELASTIC
          // R(9,9) = 1.0;
          // eigenvector corresponding to the additional damage variable
          real I1 = local.epsxx_alpha + local.epsyy_alpha + local.epszz_alpha;
          real I2 = local.epsxx_alpha*local.epsxx_alpha
            + local.epsyy_alpha*local.epsyy_alpha
            + local.epszz_alpha*local.epszz_alpha
            + 2*local.epsxy_alpha*local.epsxy_alpha
            + 2*local.epsyz_alpha*local.epsyz_alpha
            + 2*local.epszx_alpha*local.epszx_alpha;

          real xi;
          if (I2 > 1e-30){
            xi = I1 / std::sqrt(I2);
          } else{
            xi = 0.0;
          }

          real xiInv;
          if ( std::abs(xi) > 1e-1){
            xiInv = 1 / xi;
          } else{
            xiInv = 0.0;
          }

          
          R(9,9) = local.rho;
          R(10,10) = local.rho;
         #endif

         //===============Added for free surface BC of the strain-vel case=====================
         //The input of getTransposedFreeSurfaceGodunovState() is changed from R to R_sig
         Matrix1010 C = Matrix1010::Zero();

         C(0,0) = local.lambda + 2.0*local.mu; C(0,1) = local.lambda; C(0,2) = local.lambda;
         C(1,0) = local.lambda; C(1,1) = local.lambda + 2.0*local.mu; C(1,2) = local.lambda;
         C(2,0) = local.lambda; C(2,1) = local.lambda; C(2,2) = local.lambda + 2.0*local.mu;
         C(3,3) = 2.0*local.mu; C(4,4) = 2.0*local.mu; C(5,5) = 2.0*local.mu;
         C(6,6) = 1; C(7,7) = 1; C(8,8) = 1;

         #if USE_DAMAGEDELASTIC
          C(9,9) = 1.0;
          C(10,10) = 1.0;
         #endif

         Matrix1010 R_sig = (C*R).eval();


         if (faceType == FaceType::freeSurface) {
          MaterialType materialtype = testIfAcoustic(local.mu) ? MaterialType::acoustic : MaterialType::elastic;
          getTransposedFreeSurfaceGodunovState(materialtype, QgodLocal, QgodNeighbor, R_sig);


          Matrix1010 Qgod = Matrix1010::Zero();
          for (unsigned i = 0; i < Qgod.cols(); ++i) {
            for (unsigned j = 0; j < Qgod.rows(); ++j) {
              Qgod(i,j) = -QgodLocal(j,i);
              Qgod(i,j) = QgodLocal(j,i);
            }
          }
          Matrix1010 Qgod_temp = Matrix1010::Zero();
          Qgod_temp = ((C.inverse()*Qgod)*C).eval();

          for (unsigned i = 0; i < Qgod.cols(); ++i) {
            for (unsigned j = 0; j < Qgod.rows(); ++j) {
              QgodLocal(i,j) = -Qgod_temp(j,i);
              QgodLocal(i,j) = Qgod_temp(j,i);
            }
          }
         //===============Added for free surface BC of the strain-vel case=====================

         } else {
           Matrix1010 chi = Matrix1010::Zero();
           if (!testIfAcoustic(local.mu)) {
             chi(2,2) = 1.0;
             chi(1,1) = 1.0;

           }
           chi(0,0) = 1.0;

           const auto godunov = ((R*chi)*R.inverse()).eval();

           // QgodLocal = I - QgodNeighbor
           for (unsigned i = 0; i < godunov.cols(); ++i) {
             for (unsigned j = 0; j < godunov.rows(); ++j) {
               QgodLocal(i,j) = -godunov(j,i);
               QgodNeighbor(i,j) = godunov(j,i);
             }
           }
           for (unsigned idx = 0; idx < 11; ++idx) {
             QgodLocal(idx,idx) += 1.0;
           }
         }
      }
  }
}
#endif
