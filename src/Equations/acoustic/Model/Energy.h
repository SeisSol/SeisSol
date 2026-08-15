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
#include "Kernels/Precision.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<AcousticMaterial> {
  static constexpr auto Energies = AcousticEnergies;
  static constexpr std::size_t EnergyCount = Energies.size();
  static_assert(detail::descriptorsWellFormed(Energies),
                "energy descriptors must be named, unique, and grouped consistently");

  // output positions, looked up by name so that reordering cannot misplace a value
  static constexpr auto AcousticEnergyIdx = detail::indexOf(Energies, "acoustic_energy");
  static constexpr auto AcousticKineticIdx = detail::indexOf(Energies, "acoustic_kinetic_energy");
  static_assert(AcousticEnergyIdx < EnergyCount, "AcousticEnergy missing from the descriptor list");
  static_assert(AcousticKineticIdx < EnergyCount,
                "AcousticKinetic missing from the descriptor list");

  /// No anelastic variables. See the viscoelastic specialization for the
  /// non-trivial case; the argument is accepted uniformly so that
  /// EnergyOutput does not need to branch on the material.
  struct Moments {};
  static Moments computeMoments(const real* /*dofs*/, const real* /*dofsAne*/) { return {}; }

  static AcousticMaterial::EnergyData initEnergyData(const AcousticMaterial& /*material*/) {
    return {};
  }

  template <typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount>
      computeEnergies(const AcousticMaterial& material,
                      const AcousticMaterial::EnergyData& /*data*/,
                      const LinearViewT& /*linSub*/,
                      const QuadraticViewT& quadSub,
                      const Moments& /*moments*/,
                      std::size_t /*sim*/) {
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

    output[AcousticEnergyIdx] = curAcousticEnergy;
    output[AcousticKineticIdx] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_ACOUSTIC_MODEL_ENERGY_H_
