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

class NoFault : public BaseFrictionSolver<NoFault> {
  public:
  explicit NoFault(const FrictionLawParameters& drParameters)
      : BaseFrictionSolver<NoFault>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {}

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {
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
  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& __restrict ctx,
                                                     real time) {}

  SEISSOL_DEVICE static void preHook(FrictionLawContext& __restrict ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext& __restrict ctx) {}

  protected:
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
