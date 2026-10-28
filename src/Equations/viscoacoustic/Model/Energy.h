// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_ENERGY_H_

#include "Common/Constants.h"
#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/EnergyBase.h"
#include "Equations/viscoacoustic/Model/Datastructures.h"
#include "Model/Common.h"

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "Common/ConfigDispatch.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Solver/MultipleSimulations.h"
#endif

#include <array>
#include <cmath>
#include <cstddef>

namespace seissol::model {

template <std::size_t Mechanisms>
struct EnergyCompute<ViscoAcousticMaterial<Mechanisms>> {
  using ViscoMaterial = ViscoAcousticMaterial<Mechanisms>;

  static constexpr auto Energies =
      detail::concat(MomentumEnergies, AcousticEnergies, ViscoacousticEnergies, ViscoEnergies);
  static constexpr std::size_t EnergyCount = Energies.size();
  static_assert(detail::descriptorsWellFormed(Energies),
                "energy descriptors must be named, unique, and grouped consistently");

  // output positions, looked up by name so that reordering cannot misplace a value
  static constexpr auto MomentumXIdx = detail::indexOf(Energies, "momentumX");
  static constexpr auto MomentumYIdx = detail::indexOf(Energies, "momentumY");
  static constexpr auto MomentumZIdx = detail::indexOf(Energies, "momentumZ");
  static constexpr auto AcousticPotentialIdx =
      detail::indexOf(Energies, "acoustic_potential_energy");
  static constexpr auto AcousticKineticIdx = detail::indexOf(Energies, "acoustic_kinetic_energy");
  static constexpr auto AnelasticPotentialIdx =
      detail::indexOf(Energies, "anelastic_potential_energy");
  static constexpr auto ViscousDissipationIdx =
      detail::indexOf(Energies, "viscous_dissipation_rate");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(AcousticPotentialIdx < EnergyCount,
                "AcousticPotential missing from the descriptor list");
  static_assert(AcousticKineticIdx < EnergyCount,
                "AcousticKinetic missing from the descriptor list");
  static_assert(AnelasticPotentialIdx < EnergyCount,
                "AnelasticPotential missing from the descriptor list");
  static_assert(ViscousDissipationIdx < EnergyCount,
                "ViscousDissipation missing from the descriptor list");

  /**
   * Cell moments involving the anelastic variables, see codegen kernels
   * `momentQaneQaneCompute` and `momentQQaneCompute`, for a configuration with the anelastic
   * solver. With the solver linearck, the anelastic variables are among the unknowns, so that their
   * moments are second moments of the unknowns, which computeEnergies gets anyway.
   */
  template <typename Cfg, bool Anelastic = Cfg::Solver == SolverType::LinearCKAnelastic>
  struct Moments {};

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
  template <typename Cfg>
  struct Moments<Cfg, true> {
    alignas(Alignment) Real<Cfg> ane[tensor::momentQaneQane<Cfg>::size()]{};
    alignas(Alignment) Real<Cfg> cross[tensor::momentQQane<Cfg>::size()]{};
  };
#endif

  template <typename Cfg>
  static Moments<Cfg> computeMoments([[maybe_unused]] const Real<Cfg>* dofs,
                                     [[maybe_unused]] const Real<Cfg>* dofsAne) {
    Moments<Cfg> moments{};
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
    if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
      constexpr auto Variant = configIdOf<Cfg>();

      runtime::kernel::momentQaneQaneCompute aneKrnl;
      aneKrnl.Qane = runtime::init::Qane::view(Variant, dofsAne);
      aneKrnl.momentQaneQane = runtime::init::momentQaneQane::view(Variant, moments.ane);
      aneKrnl.execute(Variant);

      runtime::kernel::momentQQaneCompute crossKrnl;
      crossKrnl.Q = runtime::init::Q::view(Variant, dofs);
      crossKrnl.Qane = runtime::init::Qane::view(Variant, dofsAne);
      crossKrnl.momentQQane = runtime::init::momentQQane::view(Variant, moments.cross);
      crossKrnl.execute(Variant);
    }
#endif
    return moments;
  }

  template <typename Cfg>
  static typename ViscoMaterial::template EnergyData<Cfg>
      initEnergyData(const ViscoMaterial& /*material*/) {
    return {};
  }

  /**
   * The viscoelastic system is a generalized Maxwell body,
   *
   *   dsigma/dt        = C_u : epsdot - sum_l D^(l) : vartheta^(l)
   *   dvartheta^(l)/dt = omega_l ( epsdot - vartheta^(l) )
   *
   * with D^(l) the isotropic modulus defect of mechanism l (== -E, see
   * getTransposedSourceCoefficientTensor) and C_u the unrelaxed stiffness, i.e.
   * material.lambda/mu after fitAttenuation.
   *
   * The anelastic variables carry the dimension of a strain *rate* -- the
   * relaxation equation forces [vartheta] = [epsdot]. Writing xi^(l) for their
   * time integral and e^(l) := eps - xi^(l) for the spring strain of branch l,
   * one gets xidot^(l) = omega_l e^(l), hence
   *
   *   e^(l) = vartheta^(l) / omega_l
   *
   * so no history is needed: both the stored branch energy and the dissipation
   * rate are functions of the instantaneous state. The total strain follows from
   *
   *   sigma_r := sigma - sum_l D^(l) vartheta^(l) / omega_l  =  C_r : eps
   *
   * with C_r = C_u - sum_l D^(l) the relaxed stiffness. Then
   *
   *   W         = 1/2 eps:C_r:eps + 1/2 sum_l e^(l):D^(l):e^(l)
   *   Wdot_diss = sum_l (1/omega_l) vartheta^(l):D^(l):vartheta^(l)  >= 0
   *
   * Note that 1/2 sigma:C_u^-1:sigma -- the elastic formula -- is *not* the
   * stored energy, and neither is 1/2 sigma:C_r^-1:sigma; the sigma_r correction
   * cannot be absorbed into a choice of moduli.
   *
   * Everything below evaluates these in a volumetric/deviatoric split, which
   * avoids inverting any 6x6 matrix.
   */
  template <typename Cfg, typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount>
      computeEnergies(const ViscoMaterial& material,
                      const typename ViscoMaterial::template EnergyData<Cfg>& /*data*/,
                      const LinearViewT& linSub,
                      const QuadraticViewT& quadSub,
                      [[maybe_unused]] const Moments<Cfg>& moments,
                      [[maybe_unused]] std::size_t sim) {
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
    if constexpr (Cfg::Solver == SolverType::LinearCKAnelastic) {
      const auto aneFused = init::momentQaneQane<Cfg>::view::create(moments.ane);
      const auto crossFused = init::momentQQane<Cfg>::view::create(moments.cross);
      const auto ane = multisim::simtensor<Cfg>(aneFused, sim);
      const auto cross = multisim::simtensor<Cfg>(crossFused, sim);
      return computeEnergiesFrom(material, linSub, quadSub, cross, ane);
    } else
#endif
    {
      // the anelastic variables follow the elastic quantities among the unknowns, a block per
      // mechanism
      const auto anelasticIndex = [](std::size_t component, std::size_t mechanism) {
        return ViscoMaterial::NumElasticQuantities + mechanism * ViscoMaterial::NumberPerMechanism +
               component;
      };
      const auto cross = [&](std::size_t i, std::size_t j, std::size_t m) {
        return quadSub(i, anelasticIndex(j, m));
      };
      const auto ane = [&](std::size_t i, std::size_t j, std::size_t m, std::size_t n) {
        return quadSub(anelasticIndex(i, m), anelasticIndex(j, n));
      };
      return computeEnergiesFrom(material, linSub, quadSub, cross, ane);
    }
  }

  /**
   * The energies of a cell, given the moments of its anelastic variables: cross(i, j, m) the one of
   * quantity i and component j of mechanism m, ane(i, j, m, n) the one of component i of mechanism
   * m and component j of mechanism n.
   */
  template <typename LinearViewT, typename QuadraticViewT, typename CrossT, typename AneT>
  static std::array<double, EnergyCount> computeEnergiesFrom(const ViscoMaterial& material,
                                                             const LinearViewT& linSub,
                                                             const QuadraticViewT& quadSub,
                                                             const CrossT& cross,
                                                             const AneT& ane) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = ViscoMaterial::VelocityOffset;
    const auto rho = material.getDensity();

    const auto u = linSub(0, UIdx + 0);
    const auto v = linSub(0, UIdx + 1);
    const auto w = linSub(0, UIdx + 2);
    const auto uu = quadSub(UIdx + 0, UIdx + 0);
    const auto vv = quadSub(UIdx + 1, UIdx + 1);
    const auto ww = quadSub(UIdx + 2, UIdx + 2);
    const double curKineticEnergy = 0.5 * rho * (uu + vv + ww);

    // the momentum is material-independent; in particular it must also be reported for
    // acoustic cells (otherwise mixed acoustic/elastic setups silently lose momentum)
    output[MomentumXIdx] = rho * u;
    output[MomentumYIdx] = rho * v;
    output[MomentumZIdx] = rho * w;

    // Moments of the traces. The quantity is a single isotropic stress
    // component, so the trace of the tensor it stands for is three times over
    // it -- which is what the weights below carry, and what keeps the
    // expressions further down identical to the ones for a stress tensor.
    constexpr auto Weight = traceWeight(roleKind(ViscoMaterial::PrimaryGroups, FaceRole::Traction));
    const auto traceSigmaSq = Weight * Weight * quadSub(0, 0);
    const auto traceSigmaTheta = [&](std::size_t m) { return Weight * cross(0, 0, m); };
    // No weight here: the memory variable is driven by the divergence of the
    // velocity, so it already is the trace of its tensor counterpart rather
    // than one component of it.
    const auto traceThetaTheta = [&](std::size_t m, std::size_t n) { return ane(0, 0, m, n); };

    // sigma_r = sigma - sum_m D^(m) vartheta^(m) / omega_m, split into its
    // volumetric and deviatoric parts:
    //   tr  sigma_r = tr sigma  - sum_m 3 dK_m  tr vartheta^(m)  / omega_m
    //   dev sigma_r = dev sigma - sum_m 2 dMu_m dev vartheta^(m) / omega_m
    std::array<double, zeroGuard(Mechanisms)> volFactor{};
    for (std::size_t m = 0; m < Mechanisms; ++m) {
      const auto inverseOmega = 1.0 / material.omega[m];
      volFactor[m] = 3.0 * material.getDeltaBulk(m) * inverseOmega;
    }

    // Volumetric and deviatoric parts are orthogonal, so sigma_r splits cleanly:
    //   tr  sigma_r = tr sigma  - sum_m 3 dK_m  tr vartheta^(m)  / omega_m
    //   dev sigma_r = dev sigma - sum_m 2 dMu_m dev vartheta^(m) / omega_m
    double traceSigmaRSq = traceSigmaSq;
    for (std::size_t m = 0; m < Mechanisms; ++m) {
      traceSigmaRSq -= 2.0 * volFactor[m] * traceSigmaTheta(m);
      for (std::size_t n = 0; n < Mechanisms; ++n) {
        traceSigmaRSq += volFactor[m] * volFactor[n] * traceThetaTheta(m, n);
      }
    }

    // branch springs: e^(m) = vartheta^(m) / omega_m
    double branchEnergy = 0.0;
    double dissipationRate = 0.0;
    for (std::size_t m = 0; m < Mechanisms; ++m) {
      const auto inverseOmega = 1.0 / material.omega[m];
      const auto traceSq = traceThetaTheta(m, m) * inverseOmega * inverseOmega;
      const auto quadratic = material.getDeltaBulk(m) * traceSq;
      branchEnergy += 0.5 * quadratic;
      dissipationRate += material.omega[m] * quadratic;
    }

    output[AnelasticPotentialIdx] = branchEnergy;
    output[ViscousDissipationIdx] = dissipationRate;

    // No shear branch survives (dMu_m = mu_u * beta_m vanishes with mu_u), so
    // only the volumetric part contributes.
    const auto bulkRelaxed = material.getBulkRelaxed();
    output[AcousticPotentialIdx] = traceSigmaRSq / (18.0 * bulkRelaxed);
    output[AcousticKineticIdx] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_VISCOACOUSTIC_MODEL_ENERGY_H_
