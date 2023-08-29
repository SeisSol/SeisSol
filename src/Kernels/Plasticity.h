/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de, https://www.geophysik.uni-muenchen.de/Members/wollherr)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 * Plasticity kernel of SeisSol.
 **/

#ifndef KERNELS_PLASTICITY_H_
#define KERNELS_PLASTICITY_H_

#include "Model/plasticity.hpp"
#include <Initializer/typedefs.hpp>
#include <Initializer/BatchRecorders/DataTypes/ConditionalTable.hpp>
#include <limits>
#include "Common/configtensor.hpp"

namespace seissol::kernels {

template<typename Config>
class Plasticity {
public:
  using RealT = typename Config::RealT;
  /** Returns 1 if there was plastic yielding otherwise 0.
   */
  static unsigned computePlasticity( double                      oneMinusIntegratingFactor,
                                     double                      timeStepWidth,
                                     double                      T_v,
                                     GlobalData<Config> const*           global,
                                     seissol::model::PlasticityData<RealT> const*       plasticityData,
                                     RealT                        degreesOfFreedom[Yateto<Config>::Tensor::Q::size()],
                                     RealT*                       pstrain);

  static unsigned computePlasticityBatched(double relaxTime,
                                           double timeStepWidth,
                                           double T_v,
                                           GlobalData<Config> const *global,
                                           initializers::recording::ConditionalPointersToRealsTable &table,
                                           seissol::model::PlasticityData<RealT> *plasticity);

  static void flopsPlasticity(  long long&  o_nonZeroFlopsCheck,
                                long long&  o_hardwareFlopsCheck,
                                long long&  o_nonZeroFlopsYield,
                                long long&  o_hardwareFlopsYield );
};

}

#endif

