// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"

namespace seissol::dr::friction_law::gpu {

template <typename Cfg>
class NoFault : public BaseFrictionSolver<Cfg, NoFault<Cfg>> {
  public:
  using real = Real<Cfg>;

  // the slip rate of a sub-step is the one the friction solve finds at its time point
  static constexpr bool SlipRateAtIntervalEnds = false;

  explicit NoFault(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, NoFault<Cfg>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {}

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext<Cfg>& __restrict ctx,
                                                   uint32_t /*timeIndex*/) {
    // calculate traction
    ctx.tractionResults.traction1 = ctx.faultStresses.traction1;
    ctx.tractionResults.traction2 = ctx.faultStresses.traction2;
    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = ctx.tractionResults.traction1;
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = ctx.tractionResults.traction2;
  }

  //! no dynamic stress time
  SEISSOL_DEVICE static real dynamicStressQuantity(FrictionLawContext<Cfg>& __restrict /*ctx*/) {
    return 0;
  }
  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext<Cfg>& __restrict ctx,
                                                     real previousQuantity,
                                                     const common::SubStepTimes<Cfg>& times) {}

  SEISSOL_DEVICE static void preHook(FrictionLawContext<Cfg>& __restrict ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext<Cfg>& __restrict ctx) {}

  protected:
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
