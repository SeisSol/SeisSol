// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_

namespace seissol::dr::friction_law::cpu {
template <typename Cfg>
class NoTP {
  public:
  using real = Real<Cfg>;

  explicit NoTP(seissol::initializer::parameters::DRParameters* drParameters) {};

  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}

  void calcFluidPressure(std::array<real, misc::NumPaddedPoints<Cfg>>& normalStress,
                         real (*mu)[misc::NumPaddedPoints<Cfg>],
                         std::array<real, misc::NumPaddedPoints<Cfg>>& slipRateMagnitude,
                         real deltaT,
                         bool saveTmpInTP,
                         uint32_t timeIndex,
                         std::size_t ltsFace) {}

  [[nodiscard]] static real getFluidPressure(unsigned /*unused*/, unsigned /*unused*/) {
    return 0;
  };
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_THERMALPRESSURIZATION_NOTP_H_
