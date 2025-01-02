/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de,
 *https://www.geophysik.uni-muenchen.de/Members/wollherr)
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

#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Model/Plasticity.h"
#include "Parallel/Runtime/Stream.h"
#include "generated_code/tensor.h"
#include <limits>

namespace seissol::kernels {

class Plasticity {
  public:
  /** Returns 1 if there was plastic yielding otherwise 0.
   */
  static unsigned computePlasticity(double oneMinusIntegratingFactor,
                                    double timeStepWidth,
                                    double tV,
                                    const GlobalData* global,
                                    const seissol::model::PlasticityData* plasticityData,
                                    real degreesOfFreedom[tensor::Q::size()],
                                    real* pstrain);

  static unsigned
      computePlasticityBatched(double oneMinusIntegratingFactor,
                               double timeStepWidth,
                               double tV,
                               const GlobalData* global,
                               initializer::recording::ConditionalPointersToRealsTable& table,
                               seissol::model::PlasticityData* plasticityData,
                               seissol::parallel::runtime::StreamRuntime& runtime);

  static void flopsPlasticity(long long& nonZeroFlopsCheck,
                              long long& hardwareFlopsCheck,
                              long long& nonZeroFlopsYield,
                              long long& hardwareFlopsYield);
};

} // namespace seissol::kernels

#endif
