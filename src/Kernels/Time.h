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
};

#endif

