// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ENERGY_H_

#include "Equations/elastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<ElasticMaterial> {
  static constexpr std::size_t EnergyCount = 7;
  static inline const std::array<std::string, EnergyCount> Energies{
      "momentum-x",
      "momentum-y",
      "momentum-z",
      "acoustic-energy",
      "elastic-energy",
      "acoustic-kinetic-energy",
      "elastic-kinetic-energy",
  };

  static ElasticMaterial::EnergyData initEnergyData(const ElasticMaterial& /*material*/) {
    return {};
  }

  static std::array<double, EnergyCount>
      computeEnergies(const ElasticMaterial& material,
                      const AcousticMaterial::EnergyData& /*data*/,
                      const init::massLPR::view::type& linSub,
                      const init::massSPR::view::type& quadSub) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = ElasticMaterial::TractionQuantities;
    const auto rho = material.getDensity();

    const auto u = linSub(0, UIdx + 0);
    const auto v = linSub(0, UIdx + 1);
    const auto w = linSub(0, UIdx + 2);
    const auto uu = quadSub(UIdx + 0, UIdx + 0);
    const auto vv = quadSub(UIdx + 1, UIdx + 1);
    const auto ww = quadSub(UIdx + 2, UIdx + 2);
    const double curKineticEnergy = 0.5 * rho * (uu + vv + ww);
    const double curMomentumX = rho * u;
    const double curMomentumY = rho * v;
    const double curMomentumZ = rho * w;

    if (std::abs(material.getMuBar()) < 10e-14) {
      // Acoustic
      constexpr std::size_t PIdx = 0;
      const auto k = material.getLambdaBar();
      const auto pp = quadSub(PIdx, PIdx);
      const double curAcousticEnergy = pp / (2 * k);

      output[3] = curAcousticEnergy;
      output[5] = curKineticEnergy;
    } else {
      // Elastic
      auto getStressIndex = [](int i, int j) {
        const static auto Lookup =
            std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
        return Lookup[i][j];
      };
      output[0] = curMomentumX;
      output[1] = curMomentumY;
      output[2] = curMomentumZ;

      auto getStressPair = [&](int i1, int j1, int i2, int j2) {
        return quadSub(getStressIndex(i1, j1), getStressIndex(i2, j2));
      };

      const auto lambda = material.getLambdaBar();
      const auto mu = material.getMuBar();
      const auto factor = -lambda / (2.0 * mu * (3.0 * lambda + 2.0 * mu));
      auto computeStressStrain = [&](int i, int j) {
        double stressstrain = 0.0;
        if (i == j) {
          stressstrain += factor * (getStressPair(i, j, 0, 0) + getStressPair(i, j, 1, 1) +
                                    getStressPair(i, j, 2, 2));
        }
        stressstrain += 1.0 / (2.0 * mu) * getStressPair(i, j, i, j);
        return stressstrain;
      };
      double curElasticEnergy = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          curElasticEnergy += computeStressStrain(i, j);
        }
      }

      output[4] = 0.5 * curElasticEnergy;
      output[6] = curKineticEnergy;
    }

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ENERGY_H_
