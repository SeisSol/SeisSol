/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Time kernel of SeisSol.
 **/

#include "Kernels/TimeBase.h"
#include "Kernels/Time.h"
#include "Kernels/GravitationalFreeSurfaceBC.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "DirichletBoundary.h"
#pragma GCC diagnostic pop

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include <Kernels/common.hpp>
#include <Kernels/denseMatrixOps.hpp>

#include <cstring>
#include <cassert>
#include <stdint.h>
#include <omp.h>

#include <yateto.h>

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

seissol::kernels::TimeBase::TimeBase() {
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < CONVERGENCE_ORDER; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
      // std::cout << tensor::dQ::size(order-1) << " " << yateto::computeFamilySize<tensor::dQ>() << std::endl;
    }
  }
}

void seissol::kernels::TimeBase::checkGlobalData(GlobalData const* global, size_t alignment) {
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % alignment == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % alignment == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % alignment == 0 );
}

void seissol::kernels::Time::setHostGlobalData(GlobalData const* global) {
#ifdef USE_STP
  //Note: We could use the space time predictor for elasticity.
  //This is not tested and experimental
  for (int n = 0; n < CONVERGENCE_ORDER; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        m_krnlPrototype.kDivMTSub(d,n) = init::kDivMTSub::Values[tensor::kDivMTSub::index(d,n)];
      }
    }
    m_krnlPrototype.selectModes(n) = init::selectModes::Values[tensor::selectModes::index(n)];
  }
  m_krnlPrototype.Zinv = init::Zinv::Values;
  m_krnlPrototype.timeInt = init::timeInt::Values;
  m_krnlPrototype.wHat = init::wHat::Values;
#else //USE_STP
  checkGlobalData(global, ALIGNMENT);

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
  // m_krnlNonlVolPrototype.kDivM = global->stiffnessMatrices;

  projectDerivativeToNodalBoundaryRotated.V3mTo2nFace = global->V3mTo2nFace;

#endif //USE_STP
}

void seissol::kernels::Time::setGlobalData(const CompoundGlobalData& global) {
  setHostGlobalData(global.onHost);

#ifdef ACL_DEVICE
  assert(global.onDevice != nullptr);
  const auto deviceAlignment = device.api->getGlobMemAlignment();
  checkGlobalData(global.onDevice, deviceAlignment);
  deviceKrnlPrototype.kDivMT = global.onDevice->stiffnessMatricesTransposed;
  deviceDerivativeToNodalBoundaryRotated.V3mTo2nFace = global.onDevice->V3mTo2nFace;

#endif
}

void seissol::kernels::Time::computeAder(double i_timeStepWidth,
                                         LocalData& data,
                                         LocalTmp& tmp,
                                         real o_timeIntegrated[tensor::I::size()],
                                         real* o_timeDerivatives,
                                         double startTime,
                                         bool updateDisplacement) {

  assert(reinterpret_cast<uintptr_t>(data.dofs) % ALIGNMENT == 0 );
  assert(reinterpret_cast<uintptr_t>(o_timeIntegrated) % ALIGNMENT == 0 );
  assert(o_timeDerivatives == nullptr || reinterpret_cast<uintptr_t>(o_timeDerivatives) % ALIGNMENT == 0);

  // Only a small fraction of cells has the gravitational free surface boundary condition
  updateDisplacement &= std::any_of(std::begin(data.cellInformation.faceTypes),
                                    std::end(data.cellInformation.faceTypes),
                                    [](const FaceType f) {
                                      return f == FaceType::freeSurfaceGravity;
                                    });

#ifdef USE_STP
  //Note: We could use the space time predictor for elasticity.
  //This is not tested and experimental
  alignas(PAGESIZE_STACK) real stpRhs[tensor::spaceTimePredictor::size()];
  alignas(PAGESIZE_STACK) real stp[tensor::spaceTimePredictor::size()]{};
  kernel::spaceTimePredictor krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }
  krnl.Q = const_cast<real*>(data.dofs);
  krnl.I = o_timeIntegrated;
  krnl.timestep = i_timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;
  krnl.execute();
#else //USE_STP
  // real epsInitxx = 4.63e-4; // eps_xx0
  // real epsInityy = -1.85e-3; // eps_yy0
  // real epsInitzz = 4.63e-4; // eps_zz0
  // real epsInitxy = 1.11e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // real epsInitxx = -9.26e-4; // eps_xx0
  // real epsInityy = -9.26e-4; // eps_yy0
  // real epsInitzz = -9.26e-4; // eps_zz0
  // real epsInitxy = 1.11e-3; // eps_xx0
  // real epsInityz = -0e-1; // eps_yy0
  // real epsInitzx = -0e-1; // eps_zz0

  // tpv 5
  // real epsInitxx = 3.73854e-4; // eps_xx0
  // real epsInityy = -1.4963e-3; // eps_yy0
  // real epsInitzz = 3.73854e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xx0
  // real epsInityz = -0e-1; // eps_yy0
  // real epsInitzx = -0e-1; // eps_zz0

  // real epsInitxx = -1.8738e-4; // eps_xx0
  // real epsInityy = -1.1225e-3; // eps_yy0
  // real epsInitzz = -1.8738e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // // tpv5 28.0 deg
  // real epsInitxx = 3.7986e-4; // eps_xx0
  // real epsInityy = -1.0383e-3; // eps_yy0
  // real epsInitzz = -1.0072e-3; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // // tpv5 45.0 deg
  // real epsInitxx = -7.4861e-4; // eps_xx0
  // real epsInityy = -7.4861e-4; // eps_yy0
  // real epsInitzz = -7.4861e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // tpv5 45.0 deg, xi 0.77
  // real epsInitxx = -1.0072e-3; // eps_xx0
  // real epsInityy = -1.0383e-3; // eps_yy0
  // real epsInitzz = 3.7986e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // // tpv5 44.6 deg, Zhao's
  // real epsInitxx = -9.5732e-4; // eps_xx0
  // real epsInityy = -9.8849e-4; // eps_yy0
  // real epsInitzz = 1.8035e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // // tpv5 30.9 deg, Zhao's y-x
  // real epsInitxx = 1.8035e-4; // eps_xx0
  // real epsInityy = -9.8849e-4; // eps_yy0
  // real epsInitzz = -9.5732e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // // tpv5 35.4 deg, Zhao's y-x
  // real epsInitxx = -2.9027e-4; // eps_xx0
  // real epsInityy = -1.0539e-3; // eps_yy0
  // real epsInitzz = -2.9027e-4; // eps_zz0
  // real epsInitxy = 1.0909e-3; // eps_xy0
  // real epsInityz = -0e-1; // eps_yz0
  // real epsInitzx = -0e-1; // eps_zx0

  // // zero initial stress
  //   real epsInitxx = -1.0e-20; // eps_xx0
  //   real epsInityy = -1.0e-20; // eps_yy0
  //   real epsInitzz = -1.0e-20; // eps_zz0
  //   real epsInitxy = -0.0; // eps_xy0
  //   real epsInityz = -0e-1; // eps_yz0
  //   real epsInitzx = -0e-1; // eps_zx0

  // Benchmark
  real epsInitxx = -2.81e-4; // eps_xx0
  real epsInityy = -1.06e-3; // eps_yy0
  real epsInitzz = -2.81e-4; // eps_zz0
  real epsInitxy = 1.0909e-3; // eps_xy0
  real epsInityz = -0e-1; // eps_yz0
  real epsInitzx = -0e-1; // eps_zx0


  real const damage_para1 = data.material.local.Cd; // 1.2e-4*2;

  real const break_coeff = 1e2*damage_para1;
  real const beta_alpha = 0.05;

  const real aB0 = 7.92418e9;
  const real aB1 = -22.7919e9;
  const real aB2 = 20.3222e9;
  const real aB3 = -5.25836e9;
  const real Cg = 1e-10;
  const real m1 = 10;
  const real m2 = 1;

  kernel::damageConvertToNodal d_converToKrnl;
  #ifdef USE_DAMAGEDELASTIC
  // Compute the nodal solutions
  alignas(PAGESIZE_STACK) real solNData[tensor::QNodal::size()];
  // auto solN = tensor::QNodal::view::create(solNData);
  d_converToKrnl.v = init::v::Values;
  d_converToKrnl.QNodal = solNData;
  d_converToKrnl.Q = data.dofs;
  d_converToKrnl.execute();

  // Compute rhs of damage evolution
  alignas(PAGESIZE_STACK) real fNodalData[tensor::FNodal::size()] = {0};
  alignas(PAGESIZE_STACK) real sxxNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
  alignas(PAGESIZE_STACK) real syyNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
  alignas(PAGESIZE_STACK) real szzNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
  alignas(PAGESIZE_STACK) real sxyNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
  alignas(PAGESIZE_STACK) real syzNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};
  alignas(PAGESIZE_STACK) real szxNodal[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = {0};

  real* exxNodal = (solNData + 0 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* eyyNodal = (solNData + 1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* ezzNodal = (solNData + 2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* exyNodal = (solNData + 3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* eyzNodal = (solNData + 4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* ezxNodal = (solNData + 5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* alphaNodal = (solNData + 9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  real* breakNodal = (solNData + 10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS);
  
  real alpha_ave = 0.0;
  real break_ave = 0.0;
  real w_ave = 1.0/NUMBER_OF_ALIGNED_BASIS_FUNCTIONS;
  for (unsigned int q = 0; q<NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q){
    break_ave += breakNodal[q] * w_ave;
    alpha_ave += alphaNodal[q] * w_ave;
  }



  for (unsigned int q = 0; q<NUMBER_OF_ALIGNED_BASIS_FUNCTIONS; ++q){
    real EspI = (exxNodal[q]+epsInitxx) + (eyyNodal[q]+epsInityy) + (ezzNodal[q]+epsInitzz);
    real EspII = (exxNodal[q]+epsInitxx)*(exxNodal[q]+epsInitxx)
      + (eyyNodal[q]+epsInityy)*(eyyNodal[q]+epsInityy)
      + (ezzNodal[q]+epsInitzz)*(ezzNodal[q]+epsInitzz)
      + 2*(exyNodal[q]+epsInitxy)*(exyNodal[q]+epsInitxy)
      + 2*(eyzNodal[q]+epsInityz)*(eyzNodal[q]+epsInityz)
      + 2*(ezxNodal[q]+epsInitzx)*(ezxNodal[q]+epsInitzx);
    real xi;
    if (EspII > 1e-30){
      xi = EspI / std::sqrt(EspII);
    } else{
      xi = 0.0;
    }

    // Compute alpha_{cr}
    //(TODO: modify these calculations correctly), check for xi discrepancies
    real aCR = (3.0*xi*xi - 3.0)*data.material.local.gammaR*data.material.local.gammaR
    + 6.0*xi*data.material.local.gammaR*data.material.local.xi0*data.material.local.gammaR
    + 4.0*data.material.local.xi0*data.material.local.gammaR*data.material.local.xi0*data.material.local.gammaR;

    real bCR = - (8.0*data.material.local.mu0 + 6.0*data.material.local.lambda0) * data.material.local.xi0*data.material.local.gammaR
      - xi * (xi*xi* data.material.local.lambda0 + 6.0*data.material.local.mu0) * data.material.local.gammaR;

    real cCR = 4.0 * data.material.local.mu0 * data.material.local.mu0
      + 6.0 * data.material.local.mu0 * data.material.local.lambda0;

    real alphaCR1q = ( -bCR - std::sqrt(bCR*bCR - 4.0*aCR*cCR) )/(2.0*aCR);
    real alphaCR2q = 2.0*data.material.local.mu0
      /data.material.local.gammaR/(xi+2.0*data.material.local.xi0);

    real alphaCRq = 1.0;
    if (alphaCR1q > 0.0){
      if (alphaCR2q > 0.0){
        alphaCRq = std::min(1.0,
          std::min( alphaCR1q, alphaCR2q )
        );
      }
    }

            // damage stress -> seems different?
        real mu_eff = data.material.local.mu0 - alphaNodal[q]*data.material.local.gammaR*data.material.local.xi0
            - 0.5*alphaNodal[q]*data.material.local.gammaR*xi;
        real sxx_s = data.material.local.lambda0*EspI
                      - alphaNodal[q]*data.material.local.gammaR * std::sqrt(EspII)
                      + 2*mu_eff*(exxNodal[q]+epsInitxx);
        real syy_s = data.material.local.lambda0*EspI
                      - alphaNodal[q]*data.material.local.gammaR * std::sqrt(EspII)
                      + 2*mu_eff*(eyyNodal[q]+epsInityy);

        real szz_s = data.material.local.lambda0*EspI
                      - alphaNodal[q]*data.material.local.gammaR * std::sqrt(EspII)
                      + 2*mu_eff*(ezzNodal[q]+epsInitzz);

        real sxy_s = 2*mu_eff*(exyNodal[q]+epsInitxy); // just elastic needs to be considered here.
        real syz_s = 2*mu_eff*(eyzNodal[q]+epsInityz);
        real szx_s = 2*mu_eff*(ezxNodal[q]+epsInitzx);

        // breakage stress
        real sxx_b = (2.0*aB2 + 3.0*xi*aB3)*EspI
                      + aB1 * std::sqrt(EspII)
                      + (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(exxNodal[q]+epsInitxx); // just elastic needs to be considered here.
        real syy_b = (2.0*aB2 + 3.0*xi*aB3)*EspI
                      + aB1 * std::sqrt(EspII)
                      + (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(eyyNodal[q]+epsInityy);
        real szz_b = (2.0*aB2 + 3.0*xi*aB3)*EspI
                      + aB1 * std::sqrt(EspII)
                      + (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(ezzNodal[q]+epsInitzz);

        real sxy_b = (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(exyNodal[q]+epsInitxy);
        real syz_b = (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(eyzNodal[q]+epsInityz);
        real szx_b = (2.0*aB0 + aB1*xi - aB3*xi*xi*xi)*(ezxNodal[q]+epsInitzx);

        sxxNodal[q] = (1-breakNodal[q])*sxx_s + breakNodal[q]*sxx_b;
        syyNodal[q] = (1-breakNodal[q])*syy_s + breakNodal[q]*syy_b;
        szzNodal[q] = (1-breakNodal[q])*szz_s + breakNodal[q]*szz_b;
        sxyNodal[q] = (1-breakNodal[q])*sxy_s + breakNodal[q]*sxy_b;
        syzNodal[q] = (1-breakNodal[q])*syz_s + breakNodal[q]*syz_b;
        szxNodal[q] = (1-breakNodal[q])*szx_s + breakNodal[q]*szx_b;

    real Cplas = Cg * std::pow(breakNodal[q], m1);

    real sigma_mm = (sxxNodal[q] + syyNodal[q] + szzNodal[q])/3.0;
    real s11 = sxxNodal[q] - sigma_mm;
    real s22 = syyNodal[q] - sigma_mm;
    real s33 = szzNodal[q] - sigma_mm;
    real s12 = sxyNodal[q];
    real s23 = syzNodal[q];
    real s31 = szxNodal[q];


    //TODO, calculate the deviatoric stresses here and put those values in NodatData of the fluxes
    // TODO, get the constants m2 from parameters file
    // TODO, the right value on RHS is Cplas*std::pow(s_ij, m2)
    
    fNodalData[0*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -Cplas*std::pow(s11, m2);
    fNodalData[1*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -Cplas*std::pow(s22, m2);
    fNodalData[2*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -Cplas*std::pow(s33, m2);
    fNodalData[3*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -Cplas*std::pow(s12, m2);
    fNodalData[4*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -Cplas*std::pow(s23, m2);
    fNodalData[5*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = -Cplas*std::pow(s31, m2);
    fNodalData[11*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = Cplas*std::pow(s11, m2);
    fNodalData[12*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = Cplas*std::pow(s22, m2);
    fNodalData[13*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = Cplas*std::pow(s33, m2);
    fNodalData[14*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = Cplas*std::pow(s12, m2);
    fNodalData[15*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = Cplas*std::pow(s23, m2);
    fNodalData[16*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = Cplas*std::pow(s31, m2);

    //TODO: add the healing conditions correctly
    if (xi + data.material.local.xi0 > 0) {
      if (alpha_ave < 1.0){
        if (break_ave < 1.0){
          fNodalData[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
            (1 - breakNodal[q]) * 1.0/(std::exp( (alphaCRq - alphaNodal[q])/beta_alpha ) + 1.0) * break_coeff
              * EspII * (xi + data.material.local.xi0);
          fNodalData[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
            (1 - breakNodal[q]) * damage_para1
              * EspII * (xi + data.material.local.xi0);
        }
        else{
          fNodalData[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
          fNodalData[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
        }
      }
      else{
        fNodalData[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
        fNodalData[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
      }
    } else if (alpha_ave > 5e-1) {
      fNodalData[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        0.0 * damage_para1
          * EspII * (xi + data.material.local.xi0);
      fNodalData[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] =
        0.0 * damage_para1
          * EspII * (xi + data.material.local.xi0);
    }
    else {
      fNodalData[9*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
      fNodalData[10*NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + q] = 0.0;
    }
  }

  // Convert them back to modal space
  alignas(PAGESIZE_STACK) real dQModalData[tensor::dQModal::size()];
  kernel::damageAssignFToDQ d_assignFToDQ;
  d_assignFToDQ.vInv = init::vInv::Values;
  d_assignFToDQ.dQModal = dQModalData;
  d_assignFToDQ.FNodal = fNodalData;
  d_assignFToDQ.execute();

  // Assign the modal solutions to dQ(1)

  #endif

  alignas(PAGESIZE_STACK) real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()];
  auto* derivativesBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;

  kernel::derivative krnl = m_krnlPrototype;
  krnl.dQModal = dQModalData;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }

  // Optional source term
  set_ET(krnl, get_ptr_sourceMatrix(data.localIntegration.specific));

  krnl.dQ(0) = const_cast<real*>(data.dofs);
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }

  kernel::derivativeTaylorExpansion intKrnl;
  // intKrnl.dQModal = dQModalData;
  intKrnl.I = o_timeIntegrated;
  intKrnl.dQ(0) = data.dofs;
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }
  // powers in the taylor-series expansion
  intKrnl.power = i_timeStepWidth;
  intKrnl.execute0();

  if (updateDisplacement) {
    // First derivative if needed later in kernel
    std::copy_n(data.dofs, tensor::dQ::size(0), derivativesBuffer);
  } else if (o_timeDerivatives != nullptr) {
    // First derivative is not needed here but later
    // Hence stream it out
    streamstore(tensor::dQ::size(0), data.dofs, derivativesBuffer);
  }

  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
    krnl.execute(der);

    // update scalar for this derivative
    intKrnl.power *= i_timeStepWidth / real(der+1);
    intKrnl.execute(der);
  }


  // Do not compute it like this if at interface
  // Compute integrated displacement over time step if needed.
  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      if (data.faceDisplacements[face] != nullptr
          && data.cellInformation.faceTypes[face] == FaceType::freeSurfaceGravity) {
        bc.evaluate(
            face,
            projectDerivativeToNodalBoundaryRotated,
            data.boundaryMapping[face],
            data.faceDisplacements[face],
            tmp.nodalAvgDisplacements[face].data(),
            *this,
            derivativesBuffer,
            startTime,
            i_timeStepWidth,
            data.material,
            data.cellInformation.faceTypes[face]
        );
      }
    }
  }
#endif //USE_STP
}

void seissol::kernels::Time::computeBatchedAder(double i_timeStepWidth,
                                                LocalTmp& tmp,
                                                ConditionalPointersToRealsTable &dataTable,
                                                ConditionalMaterialTable &materialTable,
                                                double startTime,
                                                bool updateDisplacement) {
#ifdef ACL_DEVICE
  kernel::gpu_derivative derivativesKrnl = deviceKrnlPrototype;
  kernel::gpu_derivativeTaylorExpansion intKrnl;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if(dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto &entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    derivativesKrnl.numElements = numElements;
    intKrnl.numElements = numElements;

    intKrnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();

    unsigned starOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
      derivativesKrnl.star(i) = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Star))->getDeviceDataPtr());
      derivativesKrnl.extraOffset_star(i) = starOffset;
      starOffset += tensor::star::size(i);
    }

    unsigned derivativesOffset = 0;
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      derivativesKrnl.dQ(i) = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();
      intKrnl.dQ(i) = const_cast<const real **>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr());

      derivativesKrnl.extraOffset_dQ(i) = derivativesOffset;
      intKrnl.extraOffset_dQ(i) = derivativesOffset;

      derivativesOffset += tensor::dQ::size(i);
    }

    // stream dofs to the zero derivative
    device.algorithms.streamBatchedData((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr(),
                                        (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr(),
                                        tensor::Q::Size,
                                        derivativesKrnl.numElements,
                                        device.api->getDefaultStream());

    const auto maxTmpMem = yateto::getMaxTmpMemRequired(intKrnl, derivativesKrnl);
    real* tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(maxTmpMem * numElements));

    intKrnl.power = i_timeStepWidth;
    intKrnl.linearAllocator.initialize(tmpMem);
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute0();

    for (unsigned Der = 1; Der < CONVERGENCE_ORDER; ++Der) {
      derivativesKrnl.linearAllocator.initialize(tmpMem);
      derivativesKrnl.streamPtr = device.api->getDefaultStream();
      derivativesKrnl.execute(Der);

      // update scalar for this derivative
      intKrnl.power *= i_timeStepWidth / real(Der + 1);
      intKrnl.linearAllocator.initialize(tmpMem);
      intKrnl.streamPtr = device.api->getDefaultStream();
      intKrnl.execute(Der);
    }
    device.api->popStackMemory();
  }

  if (updateDisplacement) {
    auto& bc = tmp.gravitationalFreeSurfaceBc;
    for (unsigned face = 0; face < 4; ++face) {
      bc.evaluateOnDevice(face,
                          deviceDerivativeToNodalBoundaryRotated,
                          *this,
                          dataTable,
                          materialTable,
                          startTime,
                          i_timeStepWidth,
                          device);
    }
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    o_nonZeroFlops  += kernel::derivative::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivative::hardwareFlops(l_derivative);

    // update of time integrated DOFs
    o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(l_derivative);
  }

}

unsigned seissol::kernels::Time::bytesAder()
{
  unsigned reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();

  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const real*                       i_timeDerivatives,
                                              real                              o_timeIntegrated[tensor::I::size()] )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real) 1;
  real l_secondTerm = (real) 1;
  real l_factorial  = (real) 1;
  // std::cout << *i_timeDerivatives << std::endl;
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;

  // for (int i_out = 0; i_out<9; ++i_out){
  //   std::cout << i_timeDerivatives[20*i_out+0] << " ";
  // }
  // std::cout << i_timeDerivatives[20*9] << " "<< std::endl;

  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }
  // std::cout << m_derivativesOffsets[0] << std::endl;

  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(der+1);

    intKrnl.power  = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::computeBatchedIntegral(double i_expansionPoint,
                                                    double i_integrationStart,
                                                    double i_integrationEnd,
                                                    const real** i_timeDerivatives,
                                                    real ** o_timeIntegratedDofs,
                                                    unsigned numElements) {
#ifdef ACL_DEVICE
  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real deltaTLower = i_integrationStart - i_expansionPoint;
  real deltaTUpper = i_integrationEnd - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real firstTerm  = static_cast<real>(1.0);
  real secondTerm = static_cast<real>(1.0);
  real factorial  = static_cast<real>(1.0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  real* tmpMem = reinterpret_cast<real*>(device.api->getStackMemory(intKrnl.TmpMaxMemRequiredInBytes * numElements));

  intKrnl.I = o_timeIntegratedDofs;

  unsigned derivativesOffset = 0;
  for (size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives;
    intKrnl.extraOffset_dQ(i) = derivativesOffset;
    derivativesOffset += tensor::dQ::size(i);
  }

  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der) {
    firstTerm *= deltaTUpper;
    secondTerm *= deltaTLower;
    factorial *= static_cast<real>(der + 1);

    intKrnl.power = firstTerm - secondTerm;
    intKrnl.power /= factorial;
    intKrnl.linearAllocator.initialize(tmpMem);
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(der);
  }
  device.api->popStackMemory();
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::computeTaylorExpansion( real         time,
                                                     real         expansionPoint,
                                                     real const*  timeDerivatives,
                                                     real         timeEvaluated[tensor::Q::size()] ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeEvaluated)    % ALIGNMENT == 0 );

  // assert that this is a forward evaluation in time
  assert( time >= expansionPoint );

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;

  // iterate over time derivatives
  for(int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative+1);
  }
}

void seissol::kernels::Time::computeBatchedTaylorExpansion(real time,
                                                           real expansionPoint,
                                                           real** timeDerivatives,
                                                           real** timeEvaluated,
                                                           size_t numElements) {
#ifdef ACL_DEVICE
  assert( timeDerivatives != nullptr );
  assert( timeEvaluated != nullptr );
  assert( time >= expansionPoint );
  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");
  static_assert(kernel::gpu_derivativeTaylorExpansion::TmpMaxMemRequiredInBytes == 0);

  kernel::gpu_derivativeTaylorExpansion intKrnl;
  intKrnl.numElements = numElements;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = const_cast<const real **>(timeDerivatives);
    intKrnl.extraOffset_dQ(i) = m_derivativesOffsets[i];
  }

  // iterate over time derivatives
  const real deltaT = time - expansionPoint;
  intKrnl.power = 1.0;
  for(int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.streamPtr = device.api->getDefaultStream();
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / static_cast<real>(derivative + 1);
  }
#else
  assert(false && "no implementation provided");
#endif
}

void seissol::kernels::Time::computeDerivativeTaylorExpansion(real time,
                                                     real expansionPoint,
                                                     real const*  timeDerivatives,
                                                     real timeEvaluated[tensor::Q::size()],
                                                     unsigned order) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeEvaluated)    % ALIGNMENT == 0 );

  // assert that this is a forward evaluation in time
  assert( time >= expansionPoint );

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;

  // iterate over time derivatives
  for(unsigned derivative = order; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative+1);
  }
}


void seissol::kernels::Time::flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;

  // interate over derivatives
  for (unsigned der = 0; der < CONVERGENCE_ORDER; ++der) {
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}

unsigned int* seissol::kernels::Time::getDerivativesOffsets() {
  return m_derivativesOffsets;
}
