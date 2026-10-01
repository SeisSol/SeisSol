// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"

namespace seissol::dr::friction_law::cpu {
class NoTP {
  public:
  explicit NoTP(const FrictionLawParameters& drParameters) {};

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}

  void prepareFluidPressure(real deltaT, std::size_t ltsFace) {}

  void applyShearHeating(const std::array<real, misc::NumPaddedPoints>& normalStress,
                         const real (*mu)[misc::NumPaddedPoints],
                         const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
                         std::size_t ltsFace) {}

  void finalizeFluidPressure(const std::array<real, misc::NumPaddedPoints>& normalStress,
                             const real (*mu)[misc::NumPaddedPoints],
                             const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
                             real deltaT,
                             std::size_t ltsFace) {}

  [[nodiscard]] static real getFluidPressure(std::size_t /*unused*/, std::uint32_t /*unused*/) {
    return 0;
  };

  [[nodiscard]] static real fluidPressureOffset(std::uint32_t /*unused*/) { return 0; }
  [[nodiscard]] static real fluidPressureSlope(std::uint32_t /*unused*/) { return 0; }
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
