// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_ENERGY_H_

#include "Equations/Energy.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"
#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"

namespace seissol::model {

template <>
struct EnergyCompute<PoroElasticMaterial> {
  static constexpr auto Energies = detail::concat(MomentumEnergies, ElasticEnergies, DarcyEnergies);
  static constexpr std::size_t EnergyCount = Energies.size();
  static_assert(detail::descriptorsWellFormed(Energies),
                "energy descriptors must be named, unique, and grouped consistently");

  // output positions, looked up by name so that reordering cannot misplace a value
  static constexpr auto MomentumXIdx = detail::indexOf(Energies, "momentumX");
  static constexpr auto MomentumYIdx = detail::indexOf(Energies, "momentumY");
  static constexpr auto MomentumZIdx = detail::indexOf(Energies, "momentumZ");
  static constexpr auto ElasticEnergyIdx = detail::indexOf(Energies, "elastic_energy");
  static constexpr auto ElasticKineticIdx = detail::indexOf(Energies, "elastic_kinetic_energy");
  static constexpr auto DarcyDissipationIdx = detail::indexOf(Energies, "darcy_dissipation_rate");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(ElasticEnergyIdx < EnergyCount, "ElasticEnergy missing from the descriptor list");
  static_assert(ElasticKineticIdx < EnergyCount, "ElasticKinetic missing from the descriptor list");
  static_assert(DarcyDissipationIdx < EnergyCount,
                "DarcyDissipation missing from the descriptor list");
  static_assert(MomentumYIdx == MomentumXIdx + 1 && MomentumZIdx == MomentumXIdx + 2,
                "the momentum components must be contiguous for the loop below");

  /// No anelastic variables. See the viscoelastic specialization for the
  /// non-trivial case; the argument is accepted uniformly so that
  /// EnergyOutput does not need to branch on the material.
  struct Moments {};
  static Moments computeMoments(const real* /*dofs*/, const real* /*dofsAne*/) { return {}; }

  static PoroElasticMaterial::EnergyData initEnergyData(const PoroElasticMaterial& /*material*/) {
    return {};
  }

  template <typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount>
      computeEnergies(const PoroElasticMaterial& material,
                      const PoroElasticMaterial::EnergyData& /*data*/,
                      const LinearViewT& linSub,
                      const QuadraticViewT& quadSub,
                      const Moments& /*moments*/,
                      std::size_t /*sim*/) {
    std::array<double, EnergyCount> output{};

    // Quantity layout: 0..5 stress, 6..8 solid velocity v_s, 9 pressure,
    // 10..12 relative filtration (Darcy) velocity w.
    //
    // Note that quantities 10..12 are w, *not* the fluid velocity v_f, even
    // though Quantities names them "v1_f". This follows from the storage
    // equation encoded in getTransposedCoefficientMatrix:
    //   A(9,6) = M alpha,  A(9,10) = M   =>   dp/dt = -M ( alpha div v_s + div w )
    // which is Biot's dp/dt = -M (alpha div v_s + div w). Were quantity 10 the
    // fluid velocity, the two coefficients would have to be M(alpha - phi) and
    // M phi instead.
    constexpr auto SolidVelocityIdx = PoroElasticMaterial::VelocityOffset; // 6
    constexpr auto PressureIdx = SolidVelocityIdx + 3;                     // 9
    constexpr auto DarcyVelocityIdx = PressureIdx + 1;                     // 10

    const auto params = getAdditionalParameters(material);

    const auto vv = [&](int i) { return quadSub(SolidVelocityIdx + i, SolidVelocityIdx + i); };
    const auto vw = [&](int i) { return quadSub(SolidVelocityIdx + i, DarcyVelocityIdx + i); };
    const auto ww = [&](int i) { return quadSub(DarcyVelocityIdx + i, DarcyVelocityIdx + i); };

    double solidSq = 0.0;
    double crossTerm = 0.0;
    double darcySq = 0.0;
    for (int i = 0; i < 3; ++i) {
      solidSq += vv(i);
      crossTerm += vw(i);
      darcySq += ww(i);
    }

    // Biot kinetic energy density,
    //   T = 1/2 ( rhoBar |v_s|^2 + 2 rho_f v_s.w + m |w|^2 ),
    // which is exactly the form whose time derivative reproduces the two
    // momentum equations.
    const double curKineticEnergy =
        0.5 * (params.rhoBar * solidSq + 2 * material.rhoFluid * crossTerm + params.m * darcySq);

    // Darcy drag dissipates (viscosity / permeability) |w|^2. Reported as an
    // instantaneous rate; see the note on EnergyData in Equations/Energy.h.
    output[DarcyDissipationIdx] = material.viscosity / material.permeability * darcySq;

    // Total momentum of the mixture: rhoBar v_s + rho_f w.
    for (int i = 0; i < 3; ++i) {
      output[MomentumXIdx + i] = params.rhoBar * linSub(0, SolidVelocityIdx + i) +
                                 material.rhoFluid * linSub(0, DarcyVelocityIdx + i);
    }

    // Elastic
    auto getStressIndex = [](int i, int j) {
      const static auto Lookup =
          std::array<std::array<int, 3>, 3>{{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}}};
      return Lookup[i][j];
    };

    // Biot effective stress: sigma'_ij = sigma_ij + alpha_i delta_ij p, with p == Q[9].
    // Expanding the product of two effective stress components gives four terms.
    auto getStressPair = [&](int i1, int j1, int i2, int j2) {
      const auto s1 = getStressIndex(i1, j1);
      const auto s2 = getStressIndex(i2, j2);

      double sigma = quadSub(s1, s2);
      if (i1 == j1) {
        sigma += params.alpha(i1) * quadSub(PressureIdx, s2);
      }
      if (i2 == j2) {
        sigma += params.alpha(i2) * quadSub(s1, PressureIdx);
      }
      if (i1 == j1 && i2 == j2) {
        sigma += params.alpha(i1) * params.alpha(i2) * quadSub(PressureIdx, PressureIdx);
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
    curElasticEnergy += quadSub(PressureIdx, PressureIdx) / params.M;

    output[ElasticEnergyIdx] = 0.5 * curElasticEnergy;
    output[ElasticKineticIdx] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_ENERGY_H_
