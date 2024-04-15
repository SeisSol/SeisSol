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
#include <DynamicRupture/Typedefs.hpp>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/typedefs.hpp>
#include <Kernels/Interface.hpp>
#include <Kernels/TimeBase.h>
#include <Kernels/common.hpp>
#include <Model/datastructures.hpp>
#include <cassert>
#include <generated_code/tensor.h>
#include <kernel.h>
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
    const seissol::initializer::parameters::DamagedElasticParameters* m_damagedElasticParameters;
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
  void setDamagedElasticParameters(const seissol::initializer::parameters::DamagedElasticParameters* damagedElasticParameters){
    m_damagedElasticParameters = damagedElasticParameters;
  };

#ifdef USE_DAMAGEDELASTIC
  std::tuple<real, real, real> calculateEsp(
      const real* exxNodal,
      const real* eyyNodal,
      const real* ezzNodal,
      const real* exyNodal,
      const real* eyzNodal,
      const real* ezxNodal,
      unsigned int q,
      const seissol::initializer::parameters::DamagedElasticParameters* damagedElasticParameters);

  void calculateDynamicRuptureReceiverOutput(
      const real* dofsNPlus,
      const seissol::initializer::parameters::DamagedElasticParameters& damagedElasticParameters,
      const seissol::dr::ImpedancesAndEta* impAndEtaGet,
      real* dofsStressNPlus,
      const real* dofsNMinus,
      real* dofsStressNMinus);
  void stressToDofsDynamicRupture(real* dofsStress,
                                  const real* dofs);
  void computeNonLinearBaseFrictionLaw(const seissol::dr::ImpedancesAndEta* impAndEta,
                                       unsigned ltsFace,
                                       const real* qIPlus,
                                       real* qStressIPlus,
                                       const real* qIMinus,
                                       real* qStressIMinus);
  void updateNonLinearMaterialLocal(const real* Q_aveData,
                                    CellMaterialData* materialData,
                                    unsigned int l_cell,
                                    const std::vector<Element>& elements,
                                    const std::vector<Vertex>& vertices,
                                    unsigned int meshId,
                                    real* x,
                                    real* y,
                                    real* z);
  void updateNonLinearMaterialNeighbor(CellMaterialData* materialData,
                                       unsigned int l_cell,
                                       unsigned side,
                                       const real* Q_aveData);
  void computeNonLinearLocalIntegration(
      const seissol::kernels::LocalData& data,
      real (&QInterpolatedBodyNodal)[CONVERGENCE_ORDER][tensor::QNodal::size()],
      real (&FInterpolatedBody)[CONVERGENCE_ORDER][tensor::QNodal::size()],
      real (&sxxNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      real (&syyNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      real (&szzNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      real (&sxyNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      real (&syzNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      real (&szxNodal)[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS],
      real (&FluxInterpolatedBodyX)[CONVERGENCE_ORDER][tensor::QNodal::size()],
      real (&FluxInterpolatedBodyY)[CONVERGENCE_ORDER][tensor::QNodal::size()],
      real (&FluxInterpolatedBodyZ)[CONVERGENCE_ORDER][tensor::QNodal::size()]);
  std::tuple<real, real, real, real, real, real, real, real, real, real, real, real, real>
      calculateDamageAndBreakageStresses(real mu0,
                                         real alpha,
                                         real gammaR,
                                         real xi0,
                                         real xi,
                                         real lambda0,
                                         real EspI,
                                         real EspII,
                                         real stressXX,
                                         real stressYY,
                                         real stressZZ,
                                         real stressXY,
                                         real stressYZ,
                                         real stressXZ,
                                         real aB0,
                                         real aB1,
                                         real aB2,
                                         real aB3);
  void calculateStressesFromDamageAndBreakageStresses(real* qStressDofsXX,
                                                      real* qStressDofsYY,
                                                      real* qStressDofsZZ,
                                                      real* qStressDofsXY,
                                                      real* qStressDofsYZ,
                                                      real* qStressDofsXZ,
                                                      real* qStressDofsU,
                                                      real* qStressDofsV,
                                                      real* qStressDofsW,
                                                      real* qStressDofsDAM,
                                                      real* qStressDofsBRE,
                                                      const real* qIU,
                                                      const real* qIV,
                                                      const real* qIW,
                                                      const real* qIDAM,
                                                      const real* qIBRE,
                                                      real sxxS,
                                                      real sxxB,
                                                      real syyS,
                                                      real syyB,
                                                      real szzS,
                                                      real szzB,
                                                      real sxyS,
                                                      real sxyB,
                                                      real syzS,
                                                      real syzB,
                                                      real szxS,
                                                      real szxB,
                                                      unsigned int i);
  void updateNonLinearMaterial(seissol::model::DamagedElasticMaterial& material);
#endif
};
#endif