// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ENERGY_H_

#include "Equations/EnergyBase.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<ElasticMaterial> {
  static constexpr auto Energies =
      detail::concat(MomentumEnergies, AcousticEnergies, ElasticEnergies);
  static constexpr std::size_t EnergyCount = Energies.size();
  static_assert(detail::descriptorsWellFormed(Energies),
                "energy descriptors must be named, unique, and grouped consistently");

  // output positions, looked up by name so that reordering cannot misplace a value
  static constexpr auto MomentumXIdx = detail::indexOf(Energies, "momentumX");
  static constexpr auto MomentumYIdx = detail::indexOf(Energies, "momentumY");
  static constexpr auto MomentumZIdx = detail::indexOf(Energies, "momentumZ");
  static constexpr auto AcousticEnergyIdx = detail::indexOf(Energies, "acoustic_energy");
  static constexpr auto AcousticKineticIdx = detail::indexOf(Energies, "acoustic_kinetic_energy");
  static constexpr auto ElasticEnergyIdx = detail::indexOf(Energies, "elastic_energy");
  static constexpr auto ElasticKineticIdx = detail::indexOf(Energies, "elastic_kinetic_energy");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(AcousticEnergyIdx < EnergyCount, "AcousticEnergy missing from the descriptor list");
  static_assert(AcousticKineticIdx < EnergyCount,
                "AcousticKinetic missing from the descriptor list");
  static_assert(ElasticEnergyIdx < EnergyCount, "ElasticEnergy missing from the descriptor list");
  static_assert(ElasticKineticIdx < EnergyCount, "ElasticKinetic missing from the descriptor list");

  /// No anelastic variables. See the viscoelastic specialization for the
  /// non-trivial case; the argument is accepted uniformly so that
  /// EnergyOutput does not need to branch on the material.
  struct Moments {};
  static Moments computeMoments(const real* /*dofs*/, const real* /*dofsAne*/) { return {}; }

  static ElasticMaterial::EnergyData initEnergyData(const ElasticMaterial& /*material*/) {
    return {};
  }

  template <typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount>
      computeEnergies(const ElasticMaterial& material,
                      const ElasticMaterial::EnergyData& /*data*/,
                      const LinearViewT& linSub,
                      const QuadraticViewT& quadSub,
                      const Moments& /*moments*/,
                      std::size_t /*sim*/) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = ElasticMaterial::VelocityOffset;
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

    // the momentum is material-independent; in particular it must also be reported for
    // acoustic cells (otherwise mixed acoustic/elastic setups silently lose momentum)
    output[MomentumXIdx] = curMomentumX;
    output[MomentumYIdx] = curMomentumY;
    output[MomentumZIdx] = curMomentumZ;

    if (std::abs(material.getMuBar()) < 10e-14) {
      // Acoustic
      constexpr std::size_t PIdx = 0;
      const auto k = material.getLambdaBar();
      const auto pp = quadSub(PIdx, PIdx);
      const double curAcousticEnergy = pp / (2 * k);

      output[AcousticEnergyIdx] = curAcousticEnergy;
      output[AcousticKineticIdx] = curKineticEnergy;
    } else {
      // Elastic
      auto getStressIndex = [](int i, int j) {
        const static auto Lookup =
            std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
        return Lookup[i][j];
      };

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

      output[ElasticEnergyIdx] = 0.5 * curElasticEnergy;
      output[ElasticKineticIdx] = curKineticEnergy;
    }

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_MODEL_ENERGY_H_
