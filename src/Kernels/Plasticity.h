// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Stephanie Wollherr

#ifndef SEISSOL_SRC_KERNELS_PLASTICITY_H_
#define SEISSOL_SRC_KERNELS_PLASTICITY_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/Typedefs.h"
#include "Model/Plasticity.h"
#include "Parallel/Runtime/Stream.h"

#include <cmath>
#include <limits>

namespace seissol::kernels {

class Plasticity {
  public:
  static constexpr double computeRelaxTime(double tV, double timestep) {
    return (tV > 0.0) ? -std::expm1(-timestep / tV) : 1.0;
  }

  /** Returns 1 if there was plastic yielding otherwise 0.
   */
  static std::size_t computePlasticity(double oneMinusIntegratingFactor,
                                       double timeStepWidth,
                                       double tV,
                                       const GlobalData* global,
                                       const seissol::model::PlasticityData* plasticityData,
                                       real degreesOfFreedom[tensor::Q::size()],
                                       real* pstrain);

  static void computePlasticityBatched(double timeStepWidth,
                                       double tV,
                                       const GlobalData* global,
                                       recording::ConditionalPointersToRealsTable& table,
                                       seissol::model::PlasticityData* plasticityData,
                                       std::size_t* yieldCounter,
                                       unsigned* isAdjustableVector,
                                       seissol::parallel::runtime::StreamRuntime& runtime);

  static void flopsPlasticity(std::uint64_t& nonZeroFlopsCheck,
                              std::uint64_t& hardwareFlopsCheck,
                              std::uint64_t& nonZeroFlopsYield,
                              std::uint64_t& hardwareFlopsYield);
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_PLASTICITY_H_
