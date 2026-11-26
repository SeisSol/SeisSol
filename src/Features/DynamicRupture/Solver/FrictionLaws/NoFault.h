// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_NOFAULT_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_NOFAULT_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cstdint>

namespace seissol::dr::friction_law {

/**
 * No friction: the tractions stay the ones of the stresses, once for the host and the device (see
 * FrictionLawDriver; the solvers are cpu::NoFault and gpu::NoFault).
 */
template <typename Cfg>
class NoFault {
  public:
  using real = Real<Cfg>;

  // the slip rate of a sub-step is the one the friction solve finds at its time point
  static constexpr bool SlipRateAtIntervalEnds = false;

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* /*data*/,
                                             DynamicRupture::Layer& /*layerData*/,
                                             seissol::initializer::AllocationPlace /*place*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void preHook(FrictionLawContext<Cfg, Type>& /*ctx*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void postHook(FrictionLawContext<Cfg, Type>& /*ctx*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void updateFrictionAndSlip(FrictionLawContext<Cfg, Type>& ctx,
                                                       uint32_t /*timeIndex*/) {
    using Range = common::NumPoints<Cfg, Type>::Range;
    const auto ltsFace = ctx.ltsFace;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      const auto pointIndex = firstPoint(ctx) + index;
      point(ctx.tractionResults.traction1, index) = point(ctx.faultStresses.traction1, index);
      point(ctx.tractionResults.traction2, index) = point(ctx.faultStresses.traction2, index);
      ctx.data->traction1[ltsFace][pointIndex] = point(ctx.tractionResults.traction1, index);
      ctx.data->traction2[ltsFace][pointIndex] = point(ctx.tractionResults.traction2, index);
    }
  }

  //! no dynamic stress time
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static PointValues<Cfg, Type>
      dynamicStressQuantity(FrictionLawContext<Cfg, Type>& /*ctx*/) {
    return PointValues<Cfg, Type>{};
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      saveDynamicStressOutput(FrictionLawContext<Cfg, Type>& /*ctx*/,
                              const PointValues<Cfg, Type>& /*previousQuantity*/,
                              const common::SubStepTimes<Cfg>& /*times*/) {}
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_NOFAULT_H_
