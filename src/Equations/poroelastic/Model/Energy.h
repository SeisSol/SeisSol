// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_ENERGY_H_

#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"
#include "GeneratedCode/init.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<PoroElasticMaterial> {
  static constexpr std::size_t EnergyCount = 5;
  static inline const std::array<std::string, EnergyCount> Energies{
      "momentum-x",
      "momentum-y",
      "momentum-z",
      "elastic-energy",
      "elastic-kinetic-energy",
  };

  static PoroElasticMaterial::EnergyData initEnergyData(const PoroElasticMaterial& /*material*/) {
    return {};
  }

  static std::array<double, EnergyCount>
      computeEnergies(const PoroElasticMaterial& material,
                      const AcousticMaterial::EnergyData& /*data*/,
                      const init::massLPR::view::type& linSub,
                      const init::massSPR::view::type& quadSub) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = 6;
    constexpr auto VIdx = 10;

    const auto params = getAdditionalParameters(material);

    const auto rho = material.getDensity();

    const auto u1u1 = quadSub(UIdx + 0, UIdx + 0);
    const auto u2u2 = quadSub(UIdx + 1, UIdx + 1);
    const auto u3u3 = quadSub(UIdx + 2, UIdx + 2);
    const auto v1u1 = quadSub(UIdx + 0, VIdx + 0);
    const auto v2u2 = quadSub(UIdx + 1, VIdx + 1);
    const auto v3u3 = quadSub(UIdx + 2, VIdx + 2);
    const auto v1v1 = quadSub(VIdx + 0, VIdx + 0);
    const auto v2v2 = quadSub(VIdx + 1, VIdx + 1);
    const auto v3v3 = quadSub(VIdx + 2, VIdx + 2);

    // compute squared Darcy velocity integrals (i.e. \phi * (v_f - v_s) to the square).

    const auto w1u1 = material.porosity * (v1u1 - u1u1);
    const auto w2u2 = material.porosity * (v2u2 - u2u2);
    const auto w3u3 = material.porosity * (v3u3 - u3u3);

    const auto w1w1 = material.porosity * material.porosity * (u1u1 - 2 * v1u1 + v1v1);
    const auto w2w2 = material.porosity * material.porosity * (u2u2 - 2 * v2u2 + v2v2);
    const auto w3w3 = material.porosity * material.porosity * (u3u3 - 2 * v3u3 + v3v3);

    const double curKineticEnergy =
        0.5 * (params.rhoBar * (u1u1 + u2u2 + u3u3) + 2 * material.rhoFluid * (w1u1 + w2u2 + w3u3) +
               params.m * (w1w1 + w2w2 + w3w3));

    // TODO: the momentum will still need some adjustments
    const auto u = linSub(0, UIdx + 0);
    const auto v = linSub(0, UIdx + 1);
    const auto w = linSub(0, UIdx + 2);
    const double curMomentumX = rho * u;
    const double curMomentumY = rho * v;
    const double curMomentumZ = rho * w;

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
      double sigma = quadSub(getStressIndex(i1, j1), getStressIndex(i2, j2));
      // damp the stress (Biot-effective stress)
      if (i1 == j1) {
        sigma += params.alpha(i1) * quadSub(i2, 9);
      }
      if (i2 == j2) {
        sigma += params.alpha(i2) * quadSub(9, i1);
      }
      if (i1 == i2 && i1 == j2 && i1 == j1) {
        sigma += params.alpha(i1) * params.alpha(i2) * quadSub(9, 9);
      }
      return sigma;
    };

    // "default" isotropic elastic energy
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

    // extra pressure term
    curElasticEnergy += quadSub(9, 9) / params.M;

    output[3] = 0.5 * curElasticEnergy;
    output[4] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_ENERGY_H_
