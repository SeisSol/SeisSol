// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_THERMALPRESSURIZATION_NOTP_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cstdint>

namespace seissol::dr::friction_law {

/// no thermal pressurization: the fluid pressure stays zero (see ThermalPressurization)
template <typename Cfg>
class NoTP {
  public:
  using real = Real<Cfg>;

  /// without pressurization, a step has no coefficients
  template <common::RangeType Type>
  struct StepCoefficients {};

  static void copyStorageToLocal(FrictionLawData<Cfg>* /*data*/,
                                 DynamicRupture::Layer& /*layerData*/,
                                 seissol::initializer::AllocationPlace /*place*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real getFluidPressure(const FrictionLawContext<Cfg, Type>& /*ctx*/,
                                                  std::uint32_t /*pointIndex*/) {
    return 0;
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void prepareFluidPressure(FrictionLawContext<Cfg, Type>& /*ctx*/,
                                                      std::uint32_t /*timeIndex*/,
                                                      StepCoefficients<Type>& /*step*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real fluidPressureOffset(const StepCoefficients<Type>& /*step*/,
                                                     std::uint32_t /*pointIndex*/) {
    return 0;
  }
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real fluidPressureSlope(const StepCoefficients<Type>& /*step*/,
                                                    std::uint32_t /*pointIndex*/) {
    return 0;
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      applyShearHeating(FrictionLawContext<Cfg, Type>& /*ctx*/,
                        const PointValues<Cfg, Type>& /*normalStress*/,
                        const PointValues<Cfg, Type>& /*slipRateMagnitude*/,
                        const StepCoefficients<Type>& /*step*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      finalizeFluidPressure(FrictionLawContext<Cfg, Type>& /*ctx*/,
                            const PointValues<Cfg, Type>& /*normalStress*/,
                            const PointValues<Cfg, Type>& /*slipRateMagnitude*/,
                            std::uint32_t /*timeIndex*/,
                            const StepCoefficients<Type>& /*step*/) {}
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_THERMALPRESSURIZATION_NOTP_H_
