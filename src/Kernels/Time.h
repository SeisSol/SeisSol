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
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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

#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/TimeBase.h"
#include "generated-code/tensor.h"
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <limits>
#ifdef USE_STP
#include "Numerical/BasisFunction.h"
#include <memory>
#endif

namespace seissol::kernels {

class Time : public TimeBase {
  public:
  void setHostGlobalData(const GlobalData* global);
  void setGlobalData(const CompoundGlobalData& global);

  void computeAder(double timeStepWidth,
                   LocalData& data,
                   LocalTmp& tmp,
                   real timeIntegrated[tensor::I::size()],
                   real* timeDerivatives = nullptr,
                   bool updateDisplacement = false);

#ifdef USE_STP
  void executeSTP(double timeStepWidth,
                  LocalData& data,
                  real timeIntegrated[tensor::I::size()],
                  real* stp);
  void evaluateAtTime(
      std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions,
      const real* timeDerivatives,
      real timeEvaluated[tensor::Q::size()]);
  void flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops);

#endif
  void computeBatchedAder(double timeStepWidth,
                          LocalTmp& tmp,
                          ConditionalPointersToRealsTable& dataTable,
                          ConditionalMaterialTable& materialTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops);

  unsigned bytesAder();

  void computeIntegral(double expansionPoint,
                       double integrationStart,
                       double integrationEnd,
                       const real* timeDerivatives,
                       real timeIntegrated[tensor::I::size()]);

  void computeBatchedIntegral(double expansionPoint,
                              double integrationStart,
                              double integrationEnd,
                              const real** timeDerivatives,
                              real** timeIntegratedDofs,
                              unsigned numElements,
                              seissol::parallel::runtime::StreamRuntime& runtime);

  void computeTaylorExpansion(real time,
                              real expansionPoint,
                              const real* timeDerivatives,
                              real timeEvaluated[tensor::Q::size()]);

  void computeBatchedTaylorExpansion(real time,
                                     real expansionPoint,
                                     real** timeDerivatives,
                                     real** timeEvaluated,
                                     size_t numElements,
                                     seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops);

  unsigned int* getDerivativesOffsets();
};

} // namespace seissol::kernels

#endif
