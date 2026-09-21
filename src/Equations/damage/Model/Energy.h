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
#include "GeneratedCode/quantities.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"

#include <array>
#include <cstddef>

namespace seissol::model {

/**
 * What a damaged cell reports, taken of what it transports.
 *
 * The moments this receives are of the transported tensor, whose strain
 * columns are the perturbation the state carries and whose stress is formed
 * from the total strain -- perturbation and the strain the material starts
 * from. Every energy here is of the total strain, so the initial strain enters
 * each moment where the strain does: as a constant, integrated over the
 * reference element, whose volume is a sixth.
 *
 * The free energy is half the stress contracted with the total strain, exact
 * because the energy is homogeneous of degree two in the strain; the stress it
 * is contracted with is the one the transported tensor carries, which is its
 * projection into the element's basis, and a stress contracted with a strain in
 * that basis loses nothing to the projection. What the same strain would store
 * undamaged is a quadratic form of it. Their difference is what the damage has
 * released, and of that the onset part -- the damage times a quadratic form of
 * the strain -- is integrated by the trilinear mass form with the damage as the
 * weight, so it is exact and not a cell mean. The rest is the root part, and
 * where there is breakage it carries the granular branch's share as well.
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
  static constexpr auto FreeIdx = detail::indexOf(Energies, "damage_free_energy");
  static constexpr auto UndamagedIdx = detail::indexOf(Energies, "damage_undamaged_energy");
  static constexpr auto ReleasedIdx = detail::indexOf(Energies, "damage_released_energy");
  static constexpr auto OnsetIdx = detail::indexOf(Energies, "damage_released_onset_energy");
  static constexpr auto RootIdx = detail::indexOf(Energies, "damage_released_root_energy");
  static constexpr auto MeanDamageIdx = detail::indexOf(Energies, "mean_damage");
  static constexpr auto MeanBreakageIdx = detail::indexOf(Energies, "mean_breakage");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(KineticIdx < EnergyCount, "Kinetic missing from the descriptor list");
  static_assert(FreeIdx < EnergyCount, "Free missing from the descriptor list");
  static_assert(UndamagedIdx < EnergyCount, "Undamaged missing from the descriptor list");
  static_assert(ReleasedIdx < EnergyCount, "Released missing from the descriptor list");
  static_assert(OnsetIdx < EnergyCount, "Onset missing from the descriptor list");
  static_assert(RootIdx < EnergyCount, "Root missing from the descriptor list");
  static_assert(MeanDamageIdx < EnergyCount, "MeanDamage missing from the descriptor list");
  static_assert(MeanBreakageIdx < EnergyCount, "MeanBreakage missing from the descriptor list");

  /// The damage weighs the trilinear moment: it scales a quadratic form of the
  /// strain in the free energy.
  static constexpr std::size_t WeightColumn = generated::TransportInternalOffset;

  /// No anelastic variables. The argument is accepted uniformly so that
  /// EnergyOutput does not have to branch on the material.
  struct Moments {};
  static Moments computeMoments(const real* /*dofs*/, const real* /*dofsAne*/) { return {}; }

  static DamageMaterial::EnergyData initEnergyData(const DamageMaterial& /*material*/) {
    return {};
  }

  template <typename LinearViewT, typename QuadraticViewT, typename WeightedViewT>
  static std::array<double, EnergyCount> computeEnergies(const DamageMaterial& material,
                                                         const DamageMaterial::EnergyData& /*data*/,
                                                         const LinearViewT& linSub,
                                                         const QuadraticViewT& quadSub,
                                                         const WeightedViewT& weightedSub,
                                                         const Moments& /*moments*/,
                                                         std::size_t /*sim*/) {
    std::array<double, EnergyCount> output{};

    // Columns of the transported tensor, not of the state: the two differ for
    // this material, and the stress is only in the former.
    constexpr auto StrainIdx = generated::TransportStrainOffset;
    constexpr auto StressIdx = generated::TransportTractionOffset;
    constexpr auto VelocityIdx = generated::TransportVelocityOffset;
    constexpr auto AlphaIdx = generated::TransportInternalOffset;
    constexpr auto BreakageIdx = AlphaIdx + 1;

    // The reference tetrahedron, over which the moments are taken.
    constexpr double ReferenceVolume = 1.0 / 6.0;
    // A Voigt pair counts twice in a full contraction of two symmetric tensors.
    constexpr std::array<double, 6> Weight{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    const std::array<double, 6> initial{material.epsInitXX,
                                        material.epsInitYY,
                                        material.epsInitZZ,
                                        material.epsInitXY,
                                        material.epsInitYZ,
                                        material.epsInitXZ};

    const auto rho = material.getDensity();
    output[MomentumXIdx] = rho * linSub(0, VelocityIdx + 0);
    output[MomentumYIdx] = rho * linSub(0, VelocityIdx + 1);
    output[MomentumZIdx] = rho * linSub(0, VelocityIdx + 2);
    output[KineticIdx] = 0.5 * rho *
                         (quadSub(VelocityIdx + 0, VelocityIdx + 0) +
                          quadSub(VelocityIdx + 1, VelocityIdx + 1) +
                          quadSub(VelocityIdx + 2, VelocityIdx + 2));

    // 1/2 \int sigma : (eps + eps_0)
    double stressStrain = 0.0;
    for (std::size_t c = 0; c < 6; ++c) {
      stressStrain += Weight[c] * (quadSub(StressIdx + c, StrainIdx + c) +
                                   initial[c] * linSub(0, StressIdx + c));
    }
    const double free = 0.5 * stressStrain;

    // \int lambda0 / 2 I1^2 + mu0 I2 of the total strain
    double firstLinear = 0.0;
    double firstSquare = 0.0;
    for (std::size_t a = 0; a < 3; ++a) {
      firstLinear += linSub(0, StrainIdx + a);
      for (std::size_t b = 0; b < 3; ++b) {
        firstSquare += quadSub(StrainIdx + a, StrainIdx + b);
      }
    }
    const double firstInitial = initial[0] + initial[1] + initial[2];
    const double firstTotal = firstSquare + 2.0 * firstInitial * firstLinear +
                              firstInitial * firstInitial * ReferenceVolume;
    double secondTotal = 0.0;
    double secondWeighted = 0.0;
    for (std::size_t c = 0; c < 6; ++c) {
      const auto eps = StrainIdx + c;
      secondTotal +=
          Weight[c] * (quadSub(eps, eps) + 2.0 * initial[c] * linSub(0, eps) +
                       initial[c] * initial[c] * ReferenceVolume);
      // the same with the damage as a weight, which is the trilinear moment
      secondWeighted +=
          Weight[c] * (weightedSub(eps, eps) + 2.0 * initial[c] * quadSub(AlphaIdx, eps) +
                       initial[c] * initial[c] * linSub(0, AlphaIdx));
    }
    const double undamaged = 0.5 * material.lambda0 * firstTotal + material.mu0 * secondTotal;
    const double onset = material.gammaR * material.xi0 * secondWeighted;

    output[FreeIdx] = free;
    output[UndamagedIdx] = undamaged;
    output[ReleasedIdx] = undamaged - free;
    output[OnsetIdx] = onset;
    output[RootIdx] = undamaged - free - onset;

    output[MeanDamageIdx] = linSub(0, AlphaIdx);
    output[MeanBreakageIdx] = linSub(0, BreakageIdx);

    return output;
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_ENERGY_H_
