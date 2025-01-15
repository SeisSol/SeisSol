// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
class NoTP {
  public:
  NoTP(seissol::initializer::parameters::DRParameters* drParameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  void calcFluidPressure(real (*normalStress)[misc::NumPaddedPoints],
                         real (*mu)[misc::NumPaddedPoints],
                         real (*slipRateMagnitude)[misc::NumPaddedPoints],
                         real deltaT,
                         bool saveTmpInTP) {}

  struct Details {};
  Details getCurrentLayerDetails() { return Details{}; }
  static real getFluidPressure(Details, unsigned, unsigned) { return static_cast<real>(0.0); };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_NOTP_H_
