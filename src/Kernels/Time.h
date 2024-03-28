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

#ifndef TIME_H_
#define TIME_H_

#include <DynamicRupture/Misc.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/typedefs.hpp>
#include <Kernels/Interface.hpp>
#include <Kernels/TimeBase.h>
#include <Kernels/common.hpp>
#include <cassert>
#include <generated_code/tensor.h>
#include <limits>
#ifdef USE_STP
#include <Numerical_aux/BasisFunction.h>
#include <memory>
#endif

namespace seissol {
  namespace kernels {
    class Time;
  } // namespace kernels
} // namespace seissol

class seissol::kernels::Time : public TimeBase {
  protected:
    seissol::initializer::parameters::DamagedElasticParameters* m_damagedElasticParameters;
  public:
    void setHostGlobalData(GlobalData const* global);
    void setGlobalData(const CompoundGlobalData& global);

    void computeAder(double i_timeStepWidth,
                     LocalData& data,
                     LocalTmp& tmp,
                     real o_timeIntegrated[tensor::I::size()],
                     real* o_timeDerivatives = nullptr,
                     bool updateDisplacement = false);

#ifdef USE_STP
    void executeSTP( double     i_timeStepWidth,
                     LocalData& data,
                     real       o_timeIntegrated[tensor::I::size()],
                     real*      stp );
    void evaluateAtTime(std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions, real const* timeDerivatives, real timeEvaluated[tensor::Q::size()]);
    void flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops);

#endif
    void computeBatchedAder(double i_timeStepWidth,
                            LocalTmp& tmp,
                            ConditionalPointersToRealsTable &dataTable,
                            ConditionalMaterialTable &materialTable,
                            bool updateDisplacement = false);

    void flopsAder( unsigned int &o_nonZeroFlops,
                    unsigned int &o_hardwareFlops );

    unsigned bytesAder();

    void computeIntegral( double                                      i_expansionPoint,
                          double                                      i_integrationStart,
                          double                                      i_integrationEnd,
                          real const*                                 i_timeDerivatives,
                          real                                        o_timeIntegrated[tensor::I::size()] );

    void computeBatchedIntegral(double i_expansionPoint,
                                double i_integrationStart,
                                double i_integrationEnd,
                                const real** i_timeDerivatives,
                                real ** o_timeIntegratedDofs,
                                unsigned numElements);

    void computeTaylorExpansion( real         time,
                                 real         expansionPoint,
                                 real const*  timeDerivatives,
                                 real         timeEvaluated[tensor::Q::size()] );

    void computeDerivativeTaylorExpansion(real time,
                                          real expansionPoint,
                                          real const*  timeDerivatives,
                                          real timeEvaluated[tensor::Q::size()],
                                          unsigned derivativeOrder);


  void computeBatchedTaylorExpansion(real time,
                                     real expansionPoint,
                                     real** timeDerivatives,
                                     real** timeEvaluated,
                                     size_t numElements);

  void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops);

  unsigned int* getDerivativesOffsets();
  /**
   * sets the damaged elastic parameters pointer
  */
  void setDamagedElasticParameters(seissol::initializer::parameters::DamagedElasticParameters* damagedElasticParameters){
    m_damagedElasticParameters = damagedElasticParameters;
  };
  void calculateEps(const real* exxNodal, const real* eyyNodal, const real* ezzNodal, const real* exyNodal,
      const real* eyzNodal, const real* ezxNodal,const unsigned int& q, const seissol::initializer::parameters::DamagedElasticParameters& damagedElasticParameters,
      real& EspI, real& EspII, real&xi) {
    const real epsInitxx = damagedElasticParameters.epsInitxx;
    const real epsInityy = damagedElasticParameters.epsInityy;
    const real epsInitzz = damagedElasticParameters.epsInitzz;
    const real epsInitxy = damagedElasticParameters.epsInitxy;
    const real epsInityz = damagedElasticParameters.epsInitxx;
    const real epsInitzx = damagedElasticParameters.epsInitxx;

    EspI = (exxNodal[q] + epsInitxx) + (eyyNodal[q] + epsInityy) + (ezzNodal[q] + epsInitzz);
    EspII = (exxNodal[q] + epsInitxx) * (exxNodal[q] + epsInitxx) +
                 (eyyNodal[q] + epsInityy) * (eyyNodal[q] + epsInityy) +
                 (ezzNodal[q] + epsInitzz) * (ezzNodal[q] + epsInitzz) +
                 2 * (exyNodal[q] + epsInitxy) * (exyNodal[q] + epsInitxy) +
                 2 * (eyzNodal[q] + epsInityz) * (eyzNodal[q] + epsInityz) +
                 2 * (ezxNodal[q] + epsInitzx) * (ezxNodal[q] + epsInitzx);
    if (EspII > 1e-30) {
      xi = EspI / std::sqrt(EspII);
    } else {
      xi = 0.0;
    }
  }
void computeNonLinearRusanovFlux(const CellMaterialData* materialData, const unsigned int& l_cell, const unsigned int& side, const double* timeWeights,
  const real* qIPlus, const real* qIMinus, real* rusanovFluxP, const LocalIntegrationData* localIntegration) {
    using namespace seissol::dr::misc::quantity_indices;
    /// Checked that, after reshaping, it still uses the same memory address
    /// S4: Integration in time the Rusanov flux on surface quadrature nodes.
    const unsigned DAM = 9;
    const unsigned BRE = 10;

    const real lambda0P = materialData[l_cell].local.lambda0;
    const real mu0P = materialData[l_cell].local.mu0;
    const real rho0P = materialData[l_cell].local.rho;

    const real lambda0M = materialData[l_cell].neighbor[side].lambda0;
    const real mu0M = materialData[l_cell].neighbor[side].mu0;
    const real rho0M = materialData[l_cell].neighbor[side].rho;

    const real epsInitxx = m_damagedElasticParameters->epsInitxx;
    const real epsInityy = m_damagedElasticParameters->epsInityy;
    const real epsInitzz = m_damagedElasticParameters->epsInitzz;
    const real epsInitxy = m_damagedElasticParameters->epsInitxy;
    const real epsInityz = m_damagedElasticParameters->epsInityz;
    const real epsInitzx = m_damagedElasticParameters->epsInitzx;

    const real aB0 = m_damagedElasticParameters->aB0;
    const real aB1 = m_damagedElasticParameters->aB1;
    const real aB2 = m_damagedElasticParameters->aB2;
    const real aB3 = m_damagedElasticParameters->aB3;

    real lambdaMax = 1.0 * std::sqrt((lambda0P + 2 * mu0P) / rho0P);
    real sxxP, syyP, szzP, sxyP, syzP, szxP, sxxM, syyM, szzM, sxyM, syzM, szxM;

    for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
      auto weight = timeWeights[o];

      for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints; i++) {

        real EspIp, EspIIp, xip, EspIm, EspIIm, xim;

        calculateEps(&qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints],
                                  &qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints],
                                  i,
                                  *m_damagedElasticParameters,
                                  EspIp,
                                  EspIIp,
                                  xip);
        real alphap = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
        calculateEps(&qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints],
                                  &qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints],
                                  i,
                                  *m_damagedElasticParameters,
                                  EspIm,
                                  EspIIm,
                                  xim);
        real alpham = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + DAM*seissol::dr::misc::numPaddedPoints + i];
        real lambp = (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                         (lambda0P - alphap * materialData[l_cell].local.gammaR *
                                         (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) / std::sqrt(EspIIp)) +
                     qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (2.0 * aB2 + 3.0 * xip * aB3 +
                                          aB1 * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) / std::sqrt(EspIIp));
        real mup =
            (1 - qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                (mu0P -
                 alphap * materialData[l_cell].local.xi0 * materialData[l_cell].local.gammaR -
                 0.5 * alphap * materialData[l_cell].local.gammaR * xip) +
            qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (aB0 + 0.5 * xip * aB1 - 0.5 * xip * xip * xip * aB3);

        real lambm =
            (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                (lambda0M - alpham * materialData[l_cell].neighbor[side].gammaR *
                                (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints+i] + epsInitxx) / std::sqrt(EspIIm)) +
            qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (2.0 * aB2 + 3.0 * xim * aB3 +
                                  aB1 * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) / std::sqrt(EspIIm));

        real mum = (1 - qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i]) *
                       (mu0M -
                        alpham * materialData[l_cell].neighbor[side].xi0 *
                            materialData[l_cell].neighbor[side].gammaR -
                        0.5 * alpham * materialData[l_cell].neighbor[side].gammaR * xim) +
                   qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i] * (aB0 + 0.5 * xim * aB1 - 0.5 * xim * xim * xim * aB3);

        lambdaMax =
            std::min(std::sqrt((lambp + 2 * mup) / rho0P), std::sqrt((lambm + 2 * mum) / rho0M));

        // damage stress
        real mu_eff = materialData[l_cell].local.mu0 -
                      alphap * materialData[l_cell].local.gammaR * materialData[l_cell].local.xi0 -
                      0.5 * alphap * materialData[l_cell].local.gammaR * xip;
        real sxx_sp = materialData[l_cell].local.lambda0 * EspIp -
                      alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
                      2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_sp = materialData[l_cell].local.lambda0 * EspIp -
                      alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
                      2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_sp = materialData[l_cell].local.lambda0 * EspIp -
                      alphap * materialData[l_cell].local.gammaR * std::sqrt(EspIIp) +
                      2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_sp = 2 * mu_eff * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        // breakage stress
        real sxx_bp =
            (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_bp =
            (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_bp =
            (2.0 * aB2 + 3.0 * xip * aB3) * EspIp + aB1 * std::sqrt(EspIIp) +
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_bp =
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_bp =
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_bp =
            (2.0 * aB0 + aB1 * xip - aB3 * xip * xip * xip) * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        // damage stress minus
        mu_eff = materialData[l_cell].neighbor[side].mu0 -
                 alpham * materialData[l_cell].neighbor[side].gammaR *
                     materialData[l_cell].neighbor[side].xi0 -
                 0.5 * alpham * materialData[l_cell].neighbor[side].gammaR * xim;
        real sxx_sm = materialData[l_cell].neighbor[side].lambda0 * EspIm -
                      alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
                      2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_sm = materialData[l_cell].neighbor[side].lambda0 * EspIm -
                      alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
                      2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_sm = materialData[l_cell].neighbor[side].lambda0 * EspIm -
                      alpham * materialData[l_cell].neighbor[side].gammaR * std::sqrt(EspIIm) +
                      2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_sm = 2 * mu_eff * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        // breakage stress
        real sxx_bm =
            (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx);
        real syy_bm =
            (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy);
        real szz_bm =
            (2.0 * aB2 + 3.0 * xim * aB3) * EspIm + aB1 * std::sqrt(EspIIm) +
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz);

        real sxy_bm =
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy);
        real syz_bm =
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz);
        real szx_bm =
            (2.0 * aB0 + aB1 * xim - aB3 * xim * xim * xim) * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx);

        real breakp = qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i];
        real breakm = qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + BRE*seissol::dr::misc::numPaddedPoints + i];

        sxxP = (1 - breakp) * sxx_sp + breakp * sxx_bp;
        syyP = (1 - breakp) * syy_sp + breakp * syy_bp;
        szzP = (1 - breakp) * szz_sp + breakp * szz_bp;
        sxyP = (1 - breakp) * sxy_sp + breakp * sxy_bp;
        syzP = (1 - breakp) * syz_sp + breakp * syz_bp;
        szxP = (1 - breakp) * szx_sp + breakp * szx_bp;

        sxxM = (1 - breakm) * sxx_sm + breakm * sxx_bm;
        syyM = (1 - breakm) * syy_sm + breakm * syy_bm;
        szzM = (1 - breakm) * szz_sm + breakm * szz_bm;

        sxyM = (1 - breakm) * sxy_sm + breakm * sxy_bm;
        syzM = (1 - breakm) * syz_sm + breakm * syz_bm;
        szxM = (1 - breakm) * szx_sm + breakm * szx_bm;

        rusanovFluxP[XX*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XX*seissol::dr::misc::numPaddedPoints + i] + epsInitxx));

        rusanovFluxP[YY*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YY*seissol::dr::misc::numPaddedPoints + i] + epsInityy));

        rusanovFluxP[ZZ*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + ZZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzz));

        rusanovFluxP[XY*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XY*seissol::dr::misc::numPaddedPoints + i] + epsInitxy));

        rusanovFluxP[YZ*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ( (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + YZ*seissol::dr::misc::numPaddedPoints + i] + epsInityz));

        rusanovFluxP[XZ*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-0.5 * qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) + 0.5 * (-0.5 * qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i])) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx) -
                      0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + XZ*seissol::dr::misc::numPaddedPoints + i] + epsInitzx));

        rusanovFluxP[U*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-sxxP / rho0P) + 0.5 * (-sxxM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-sxyP / rho0P) + 0.5 * (-sxyM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-szxP / rho0P) + 0.5 * (-szxM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]) - 0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + U*seissol::dr::misc::numPaddedPoints + i]));

        rusanovFluxP[V*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-sxyP / rho0P) + 0.5 * (-sxyM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-syyP / rho0P) + 0.5 * (-syyM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-syzP / rho0P) + 0.5 * (-syzM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]) - 0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + V*seissol::dr::misc::numPaddedPoints + i]));

        rusanovFluxP[W*seissol::dr::misc::numPaddedPoints + i] +=
            weight * ((0.5 * (-szxP / rho0P) + 0.5 * (-szxM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][0] +
                      (0.5 * (-syzP / rho0P) + 0.5 * (-syzM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][1] +
                      (0.5 * (-szzP / rho0P) + 0.5 * (-szzM / rho0M)) *
                          localIntegration[l_cell].surfaceNormal[side][2] +
                      0.5 * lambdaMax * (qIPlus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]) - 0.5 * lambdaMax * (qIMinus[o*seissol::dr::misc::numQuantities*seissol::dr::misc::numPaddedPoints + W*seissol::dr::misc::numPaddedPoints + i]));
      }
    }
  }
};

#endif