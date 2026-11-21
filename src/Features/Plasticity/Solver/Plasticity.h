// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Stephanie Wollherr

#ifndef SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITY_H_
#define SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITY_H_

#include "Core/Common/Real.h"
#include "Core/Kernels/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "GeneratedCode/tensor.h"

#include <cmath>
#include <limits>

namespace seissol::kernels {

/// The plastic correction of the cells of the configuration `Cfg`.
template <typename Cfg>
class Plasticity {
  public:
  using real = Real<Cfg>;

  static constexpr double computeRelaxTime(double tV, double timestep) {
    return (tV > 0.0) ? -std::expm1(-timestep / tV) : 1.0;
  }

  /** Returns 1 if there was plastic yielding otherwise 0.
   */
  static std::size_t computePlasticity(real oneMinusIntegratingFactor,
                                       real timeStepWidth,
                                       real tV,
                                       const GlobalData<Cfg>* global,
                                       const seissol::model::PlasticityData<Cfg>* plasticityData,
                                       real degreesOfFreedom[tensor::Q<Cfg>::size()],
                                       real* pstrain);

  static void computePlasticityBatched(real timeStepWidth,
                                       real tV,
                                       const GlobalData<Cfg>* global,
                                       recording::ConditionalTable& table,
                                       seissol::model::PlasticityData<Cfg>* plasticityData,
                                       std::size_t* yieldCounter,
                                       unsigned* isAdjustableVector,
                                       seissol::parallel::runtime::StreamRuntime& runtime);

  // check, yield
  static std::pair<PerformanceEstimate, PerformanceEstimate> metrics();
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_FEATURES_PLASTICITY_SOLVER_PLASTICITY_H_
