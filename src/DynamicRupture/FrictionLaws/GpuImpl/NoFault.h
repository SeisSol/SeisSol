// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_

#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

namespace seissol::dr::friction_law::gpu {

class NoFault : public BaseFrictionSolver<NoFault> {
  public:
  NoFault(seissol::initializer::parameters::DRParameters* drParameters)
      : BaseFrictionSolver<NoFault>(drParameters) {}

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {}

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& ctx, unsigned timeIndex) {
    // calculate traction
    ctx.tractionResults.traction1[timeIndex] = ctx.faultStresses.traction1[timeIndex];
    ctx.tractionResults.traction2[timeIndex] = ctx.faultStresses.traction2[timeIndex];
    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = ctx.tractionResults.traction1[timeIndex];
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = ctx.tractionResults.traction2[timeIndex];
  }

  /*
   * output time when shear stress is equal to the dynamic stress after rupture arrived
   * currently only for linear slip weakening
   */
  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& ctx) {}

  SEISSOL_DEVICE static void preHook(FrictionLawContext& ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext& ctx) {}

  protected:
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
