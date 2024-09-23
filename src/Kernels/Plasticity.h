// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de,
 *https://www.geophysik.uni-muenchen.de/Members/wollherr)
 *
 */

#ifndef SEISSOL_SRC_KERNELS_PLASTICITY_H_
#define SEISSOL_SRC_KERNELS_PLASTICITY_H_

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
      computePlasticityBatched(double relaxTime,
                               double timeStepWidth,
                               double tV,
                               const GlobalData* global,
                               initializer::recording::ConditionalPointersToRealsTable& table,
                               seissol::model::PlasticityData* plasticity,
                               seissol::parallel::runtime::StreamRuntime& runtime);

  static void flopsPlasticity(long long& nonZeroFlopsCheck,
                              long long& hardwareFlopsCheck,
                              long long& nonZeroFlopsYield,
                              long long& hardwareFlopsYield);
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_PLASTICITY_H_
