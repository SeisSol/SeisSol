// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Stephanie Wollherr

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

#endif // SEISSOL_SRC_KERNELS_PLASTICITY_H_
