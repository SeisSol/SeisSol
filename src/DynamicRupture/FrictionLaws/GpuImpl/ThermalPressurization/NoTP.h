// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"

namespace seissol::dr::friction_law::gpu {
class NoTP {
  public:
  static void copyStorageToLocal(FrictionLawData* data, DynamicRupture::Layer& layerData) {}

  SEISSOL_DEVICE static void prepareFluidPressure(FrictionLawContext& __restrict ctx,
                                                  uint32_t timeIndex) {}

  SEISSOL_DEVICE static void applyShearHeating(FrictionLawContext& __restrict ctx) {}

  SEISSOL_DEVICE static void finalizeFluidPressure(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {}

  SEISSOL_DEVICE static real getFluidPressure(FrictionLawContext& /*unused*/) {
    return static_cast<real>(0.0);
  };

  SEISSOL_DEVICE static real fluidPressureOffset(FrictionLawContext& /*unused*/) {
    return static_cast<real>(0.0);
  }
  SEISSOL_DEVICE static real fluidPressureSlope(FrictionLawContext& /*unused*/) {
    return static_cast<real>(0.0);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_
