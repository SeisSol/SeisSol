// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

namespace seissol::dr::friction_law::gpu {
class NoTP {
  public:
  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {}

  SEISSOL_DEVICE static void
      calcFluidPressure(FrictionLawContext& ctx, int timeIndex, bool saveTmpInTP) {}

  SEISSOL_DEVICE static real getFluidPressure(FrictionLawContext& /*unused*/) {
    return static_cast<real>(0.0);
  };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_
