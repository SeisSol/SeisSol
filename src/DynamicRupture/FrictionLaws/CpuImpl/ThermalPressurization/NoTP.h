// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_

#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"

namespace seissol::dr::friction_law::cpu {
template <typename Cfg>
class NoTP {
  public:
  using real = Real<Cfg>;

  explicit NoTP(const FrictionLawParameters<Real<Cfg>>& drParameters) {};

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}

  void prepareFluidPressure(real deltaT, std::size_t ltsFace) {}

  void applyShearHeating(const std::array<real, misc::NumPaddedPoints<Cfg>>& normalStress,
                         const real (*mu)[misc::NumPaddedPoints<Cfg>],
                         const std::array<real, misc::NumPaddedPoints<Cfg>>& slipRateMagnitude,
                         std::size_t ltsFace) {}

  void finalizeFluidPressure(const std::array<real, misc::NumPaddedPoints<Cfg>>& normalStress,
                             const real (*mu)[misc::NumPaddedPoints<Cfg>],
                             const std::array<real, misc::NumPaddedPoints<Cfg>>& slipRateMagnitude,
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
