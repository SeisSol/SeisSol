// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ENERGY_H_

#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/anisotropic/Model/IntegrationData.h"
#include "GeneratedCode/init.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<AnisotropicMaterial> {
  static constexpr std::size_t EnergyCount = 5;
  static inline const std::array<std::string, EnergyCount> Energies{
      "momentum-x",
      "momentum-y",
      "momentum-z",
      "elastic-energy",
      "elastic-kinetic-energy",
  };

  static AnisotropicMaterial::EnergyData initEnergyData(const AnisotropicMaterial& material) {
    Eigen::Matrix<double, 6, 6> voigt{};

    voigt(0, 0) = material.c11;
    voigt(1, 0) = material.c12;
    voigt(2, 0) = material.c13;
    voigt(3, 0) = material.c14;
    voigt(4, 0) = material.c15;
    voigt(5, 0) = material.c16;

    voigt(0, 1) = material.c12;
    voigt(1, 1) = material.c22;
    voigt(2, 1) = material.c23;
    voigt(3, 1) = material.c24;
    voigt(4, 1) = material.c25;
    voigt(5, 1) = material.c26;

    voigt(0, 2) = material.c13;
    voigt(1, 2) = material.c23;
    voigt(2, 2) = material.c33;
    voigt(3, 2) = material.c34;
    voigt(4, 2) = material.c35;
    voigt(5, 2) = material.c36;

    voigt(0, 3) = material.c14;
    voigt(1, 3) = material.c24;
    voigt(2, 3) = material.c34;
    voigt(3, 3) = material.c44;
    voigt(4, 3) = material.c45;
    voigt(5, 3) = material.c46;

    voigt(0, 4) = material.c15;
    voigt(1, 4) = material.c25;
    voigt(2, 4) = material.c35;
    voigt(3, 4) = material.c45;
    voigt(4, 4) = material.c55;
    voigt(5, 4) = material.c56;

    voigt(0, 5) = material.c16;
    voigt(1, 5) = material.c26;
    voigt(2, 5) = material.c36;
    voigt(3, 5) = material.c46;
    voigt(4, 5) = material.c56;
    voigt(5, 5) = material.c66;

    const auto inverse = voigt.inverse().eval();

    AnisotropicEnergyData data{};

    std::copy_n(inverse.data(), data.matS.size(), data.matS.begin());

    return data;
  }

  static std::array<double, EnergyCount>
      computeEnergies(const AnisotropicMaterial& material,
                      const AnisotropicMaterial::EnergyData& data,
                      const init::massLPR::view::type& linSub,
                      const init::massSPR::view::type& quadSub) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = AnisotropicMaterial::TractionQuantities;
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

    output[0] = curMomentumX;
    output[1] = curMomentumY;
    output[2] = curMomentumZ;

    double curElasticEnergy = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        curElasticEnergy += data.matS[i + 6 * j] * quadSub(i, j);
      }
    }

    output[3] = 0.5 * curElasticEnergy;
    output[4] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ENERGY_H_
