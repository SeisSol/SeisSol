// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_THERMALPRESSURIZATION_NOTP_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_THERMALPRESSURIZATION_NOTP_H_

namespace seissol::dr::friction_law::cpu {
class NoTP {
  public:
  NoTP(seissol::initializer::parameters::DRParameters* drParameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  void calcFluidPressure(std::array<real, misc::NumPaddedPoints>& normalStress,
                         real (*mu)[misc::NumPaddedPoints],
                         std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
                         real deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace) {}

  [[nodiscard]] static real getFluidPressure(unsigned /*unused*/, unsigned /*unused*/) {
    return 0;
  };
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_THERMALPRESSURIZATION_NOTP_H_
