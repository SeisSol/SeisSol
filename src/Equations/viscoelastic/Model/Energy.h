// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_ENERGY_H_

#include "Common/Constants.h"
#include "Equations/EnergyBase.h"
#include "Equations/viscoelastic/Model/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace seissol::model {

template <std::size_t Mechanisms>
struct EnergyCompute<ViscoElasticMaterial<Mechanisms>> {
  using ViscoMaterial = ViscoElasticMaterial<Mechanisms>;

  static constexpr auto Energies =
      detail::concat(MomentumEnergies, AcousticEnergies, ElasticEnergies, ViscoelasticEnergies);
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
  static constexpr auto ViscoelasticEnergyIdx = detail::indexOf(Energies, "viscoelastic_energy");
  static constexpr auto ViscousDissipationIdx =
      detail::indexOf(Energies, "viscous_dissipation_rate");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(AcousticEnergyIdx < EnergyCount, "AcousticEnergy missing from the descriptor list");
  static_assert(AcousticKineticIdx < EnergyCount,
                "AcousticKinetic missing from the descriptor list");
  static_assert(ElasticEnergyIdx < EnergyCount, "ElasticEnergy missing from the descriptor list");
  static_assert(ElasticKineticIdx < EnergyCount, "ElasticKinetic missing from the descriptor list");
  static_assert(ViscoelasticEnergyIdx < EnergyCount,
                "ViscoelasticEnergy missing from the descriptor list");
  static_assert(ViscousDissipationIdx < EnergyCount,
                "ViscousDissipation missing from the descriptor list");

  /**
   * Cell moments involving the anelastic variables, see codegen kernels
   * `momentQaneQaneCompute` and `momentQQaneCompute`.
   */
  struct Moments {
    alignas(Alignment) real ane[tensor::momentQaneQane::size()]{};
    alignas(Alignment) real cross[tensor::momentQQane::size()]{};
  };

  static Moments computeMoments(const real* dofs, const real* dofsAne) {
    Moments moments{};

    kernel::momentQaneQaneCompute aneKrnl;
    aneKrnl.M3 = init::M3::Values;
    aneKrnl.Qane = dofsAne;
    aneKrnl.momentQaneQane = moments.ane;
    aneKrnl.execute();

    kernel::momentQQaneCompute crossKrnl;
    crossKrnl.M3 = init::M3::Values;
    crossKrnl.Q = dofs;
    crossKrnl.Qane = dofsAne;
    crossKrnl.momentQQane = moments.cross;
    crossKrnl.execute();

    return moments;
  }

  static typename ViscoMaterial::EnergyData initEnergyData(const ViscoMaterial& /*material*/) {
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
  template <typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount>
      computeEnergies(const ViscoMaterial& material,
                      const typename ViscoMaterial::EnergyData& /*data*/,
                      const LinearViewT& linSub,
                      const QuadraticViewT& quadSub,
                      const Moments& moments,
                      std::size_t sim) {
    std::array<double, EnergyCount> output{};

    constexpr auto UIdx = ViscoMaterial::TractionQuantities;
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

    const auto aneFused = init::momentQaneQane::view::create(moments.ane);
    const auto crossFused = init::momentQQane::view::create(moments.cross);
    const auto ane = multisim::simtensor(aneFused, sim);
    const auto cross = multisim::simtensor(crossFused, sim);

    // Voigt order (xx, yy, zz, xy, yz, xz). The anelastic variables carry tensor
    // components, so a double contraction weights the off-diagonals by two.
    constexpr std::array<double, 6> ContractionWeight{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};

    // moments of the traces
    const auto traceSigmaSq = [&]() {
      double sum = 0.0;
      for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          sum += quadSub(i, j);
        }
      }
      return sum;
    }();
    const auto traceSigmaTheta = [&](std::size_t m) {
      double sum = 0.0;
      for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          sum += cross(i, j, m);
        }
      }
      return sum;
    };
    const auto traceThetaTheta = [&](std::size_t m, std::size_t n) {
      double sum = 0.0;
      for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          sum += ane(i, j, m, n);
        }
      }
      return sum;
    };

    // moments of the full double contractions
    const auto contractSigmaSq = [&]() {
      double sum = 0.0;
      for (std::size_t i = 0; i < 6; ++i) {
        sum += ContractionWeight[i] * quadSub(i, i);
      }
      return sum;
    }();
    const auto contractSigmaTheta = [&](std::size_t m) {
      double sum = 0.0;
      for (std::size_t i = 0; i < 6; ++i) {
        sum += ContractionWeight[i] * cross(i, i, m);
      }
      return sum;
    };
    const auto contractThetaTheta = [&](std::size_t m, std::size_t n) {
      double sum = 0.0;
      for (std::size_t i = 0; i < 6; ++i) {
        sum += ContractionWeight[i] * ane(i, i, m, n);
      }
      return sum;
    };

    const bool acoustic = std::abs(material.getMuBar()) < 10e-14;

    // sigma_r = sigma - sum_m D^(m) vartheta^(m) / omega_m, split into its
    // volumetric and deviatoric parts:
    //   tr  sigma_r = tr sigma  - sum_m 3 dK_m  tr vartheta^(m)  / omega_m
    //   dev sigma_r = dev sigma - sum_m 2 dMu_m dev vartheta^(m) / omega_m
    std::array<double, zeroGuard(Mechanisms)> volFactor{};
    std::array<double, zeroGuard(Mechanisms)> devFactor{};
    for (std::size_t m = 0; m < Mechanisms; ++m) {
      const auto inverseOmega = 1.0 / material.omega[m];
      volFactor[m] = 3.0 * material.getDeltaBulk(m) * inverseOmega;
      devFactor[m] = 2.0 * material.getDeltaMu(m) * inverseOmega;
    }

    // Deviatoric moments follow from dev a : dev b = a:b - 1/3 tr a tr b.
    const auto devSigmaSq = contractSigmaSq - traceSigmaSq / 3.0;
    const auto devSigmaTheta = [&](std::size_t m) {
      return contractSigmaTheta(m) - traceSigmaTheta(m) / 3.0;
    };
    const auto devThetaTheta = [&](std::size_t m, std::size_t n) {
      return contractThetaTheta(m, n) - traceThetaTheta(m, n) / 3.0;
    };

    // Volumetric and deviatoric parts are orthogonal, so sigma_r splits cleanly:
    //   tr  sigma_r = tr sigma  - sum_m 3 dK_m  tr vartheta^(m)  / omega_m
    //   dev sigma_r = dev sigma - sum_m 2 dMu_m dev vartheta^(m) / omega_m
    double traceSigmaRSq = traceSigmaSq;
    double devSigmaRSq = devSigmaSq;
    for (std::size_t m = 0; m < Mechanisms; ++m) {
      traceSigmaRSq -= 2.0 * volFactor[m] * traceSigmaTheta(m);
      devSigmaRSq -= 2.0 * devFactor[m] * devSigmaTheta(m);
      for (std::size_t n = 0; n < Mechanisms; ++n) {
        traceSigmaRSq += volFactor[m] * volFactor[n] * traceThetaTheta(m, n);
        devSigmaRSq += devFactor[m] * devFactor[n] * devThetaTheta(m, n);
      }
    }

    // branch springs: e^(m) = vartheta^(m) / omega_m
    double branchEnergy = 0.0;
    double dissipationRate = 0.0;
    for (std::size_t m = 0; m < Mechanisms; ++m) {
      const auto inverseOmega = 1.0 / material.omega[m];
      const auto traceSq = traceThetaTheta(m, m) * inverseOmega * inverseOmega;
      const auto devSq = (contractThetaTheta(m, m) * inverseOmega * inverseOmega) - traceSq / 3.0;
      const auto quadratic = material.getDeltaBulk(m) * traceSq +
                             (acoustic ? 0.0 : 2.0 * material.getDeltaMu(m) * devSq);
      branchEnergy += 0.5 * quadratic;
      dissipationRate += material.omega[m] * quadratic;
    }

    output[ViscoelasticEnergyIdx] = branchEnergy;
    output[ViscousDissipationIdx] = dissipationRate;

    if (acoustic) {
      // No shear branch survives (dMu_m = mu_u * beta_m vanishes with mu_u), so
      // only the volumetric part contributes.
      const auto bulkRelaxed = material.getBulkRelaxed();
      output[AcousticEnergyIdx] = traceSigmaRSq / (18.0 * bulkRelaxed);
      output[AcousticKineticIdx] = curKineticEnergy;
    } else {
      const auto bulkRelaxed = material.getBulkRelaxed();
      const auto muRelaxed = material.getMuRelaxed();
      output[ElasticEnergyIdx] =
          traceSigmaRSq / (18.0 * bulkRelaxed) + devSigmaRSq / (4.0 * muRelaxed);
      output[ElasticKineticIdx] = curKineticEnergy;
    }

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_ENERGY_H_
