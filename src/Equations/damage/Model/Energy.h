// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_ENERGY_H_

#include "Equations/EnergyBase.h"
#include "Equations/damage/Model/Datastructures.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"

#include <array>
#include <cstddef>

namespace seissol::model {

/**
 * What a damaged cell reports.
 *
 * The momenta and the kinetic energy follow from the moments of the state
 * exactly, and so do the means of the two internal variables -- which are the
 * two numbers a run of this material is watched by.
 *
 * The stored energy does not. The free energy of the rheology carries a square
 * root of the second strain invariant, and its terms in the damage are of
 * third order in the state; the moments here reach second order in
 * polynomials. It therefore needs a nodal evaluation, of the kind the step
 * kernel already performs, and the value would reach this function through the
 * material's EnergyData rather than through the moments. What it does not have
 * yet is its expression: the coefficients of the granular branch define a
 * stress, and which potential that stress is the gradient of is a question for
 * the reference and not for this header.
 */
template <>
struct EnergyCompute<DamageMaterial> {
  static constexpr auto Energies = detail::concat(MomentumEnergies, DamageEnergies);
  static constexpr std::size_t EnergyCount = Energies.size();
  static_assert(detail::descriptorsWellFormed(Energies),
                "energy descriptors must be named, unique, and grouped consistently");

  static constexpr auto MomentumXIdx = detail::indexOf(Energies, "momentumX");
  static constexpr auto MomentumYIdx = detail::indexOf(Energies, "momentumY");
  static constexpr auto MomentumZIdx = detail::indexOf(Energies, "momentumZ");
  static constexpr auto KineticIdx = detail::indexOf(Energies, "damage_kinetic_energy");
  static constexpr auto MeanDamageIdx = detail::indexOf(Energies, "mean_damage");
  static constexpr auto MeanBreakageIdx = detail::indexOf(Energies, "mean_breakage");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(KineticIdx < EnergyCount, "Kinetic missing from the descriptor list");
  static_assert(MeanDamageIdx < EnergyCount, "MeanDamage missing from the descriptor list");
  static_assert(MeanBreakageIdx < EnergyCount, "MeanBreakage missing from the descriptor list");

  /// No anelastic variables. The argument is accepted uniformly so that
  /// EnergyOutput does not have to branch on the material.
  struct Moments {};
  static Moments computeMoments(const real* /*dofs*/, const real* /*dofsAne*/) { return {}; }

  static DamageMaterial::EnergyData initEnergyData(const DamageMaterial& /*material*/) {
    return {};
  }

  template <typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount> computeEnergies(const DamageMaterial& material,
                                                         const DamageMaterial::EnergyData& /*data*/,
                                                         const LinearViewT& linSub,
                                                         const QuadraticViewT& quadSub,
                                                         const Moments& /*moments*/,
                                                         std::size_t /*sim*/) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = DamageMaterial::VelocityOffset;
    constexpr auto AlphaIdx = DamageMaterial::VelocityOffset + 3;
    constexpr auto BreakageIdx = AlphaIdx + 1;

    const auto rho = material.getDensity();
    output[MomentumXIdx] = rho * linSub(0, UIdx + 0);
    output[MomentumYIdx] = rho * linSub(0, UIdx + 1);
    output[MomentumZIdx] = rho * linSub(0, UIdx + 2);
    output[KineticIdx] =
        0.5 * rho *
        (quadSub(UIdx + 0, UIdx + 0) + quadSub(UIdx + 1, UIdx + 1) + quadSub(UIdx + 2, UIdx + 2));

    output[MeanDamageIdx] = linSub(0, AlphaIdx);
    output[MeanBreakageIdx] = linSub(0, BreakageIdx);

    return output;
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_ENERGY_H_
