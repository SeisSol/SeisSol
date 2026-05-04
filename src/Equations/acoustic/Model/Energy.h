// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ENERGY_H_

#include "Equations/acoustic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<AcousticMaterial> {
  static constexpr std::size_t EnergyCount = 7;
  static inline const std::array<std::string, EnergyCount> Energies{
      "acoustic-energy",
      "acoustic-kinetic-energy",
  };

  static AcousticMaterial::EnergyData initEnergyData(const AcousticMaterial& /*material*/) {
    return {};
  }

  static std::array<double, EnergyCount>
      computeEnergies(const AcousticMaterial& material,
                      const AcousticMaterial::EnergyData& /*data*/,
                      const init::massLPR::view::type& linSub,
                      const init::massSPR::view::type& quadSub) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = AcousticMaterial::TractionQuantities;
    const auto rho = material.getDensity();

    const auto uu = quadSub(UIdx + 0, UIdx + 0);
    const auto vv = quadSub(UIdx + 1, UIdx + 1);
    const auto ww = quadSub(UIdx + 2, UIdx + 2);
    const double curKineticEnergy = 0.5 * rho * (uu + vv + ww);

    // Acoustic
    constexpr std::size_t PIdx = 0;
    const auto k = material.getLambdaBar();
    const auto pp = quadSub(PIdx, PIdx);
    const double curAcousticEnergy = pp / (2 * k);

    output[0] = curAcousticEnergy;
    output[1] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ENERGY_H_
