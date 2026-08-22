// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ENERGY_H_

#include "Equations/EnergyBase.h"
#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/anisotropic/Model/IntegrationData.h"
#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <utils/logger.h>

namespace seissol::model {

template <>
struct EnergyCompute<AnisotropicMaterial> {
  static constexpr auto Energies = detail::concat(MomentumEnergies, ElasticEnergies);
  static constexpr std::size_t EnergyCount = Energies.size();
  static_assert(detail::descriptorsWellFormed(Energies),
                "energy descriptors must be named, unique, and grouped consistently");

  // output positions, looked up by name so that reordering cannot misplace a value
  static constexpr auto MomentumXIdx = detail::indexOf(Energies, "momentumX");
  static constexpr auto MomentumYIdx = detail::indexOf(Energies, "momentumY");
  static constexpr auto MomentumZIdx = detail::indexOf(Energies, "momentumZ");
  static constexpr auto ElasticEnergyIdx = detail::indexOf(Energies, "elastic_energy");
  static constexpr auto ElasticKineticIdx = detail::indexOf(Energies, "elastic_kinetic_energy");
  static_assert(MomentumXIdx < EnergyCount, "MomentumX missing from the descriptor list");
  static_assert(MomentumYIdx < EnergyCount, "MomentumY missing from the descriptor list");
  static_assert(MomentumZIdx < EnergyCount, "MomentumZ missing from the descriptor list");
  static_assert(ElasticEnergyIdx < EnergyCount, "ElasticEnergy missing from the descriptor list");
  static_assert(ElasticKineticIdx < EnergyCount, "ElasticKinetic missing from the descriptor list");

  /// No anelastic variables. See the viscoelastic specialization for the
  /// non-trivial case; the argument is accepted uniformly so that
  /// EnergyOutput does not need to branch on the material.
  struct Moments {};
  static Moments computeMoments(const real* /*dofs*/, const real* /*dofsAne*/) { return {}; }

  static AnisotropicMaterial::EnergyData initEnergyData(const AnisotropicMaterial& material) {
    // The c_IJ follow *standard* Voigt numbering (1=xx, 2=yy, 3=zz, 4=yz, 5=xz,
    // 6=xy), as can be read off getTransposedCoefficientMatrix in
    // AnisotropicSetup.h. The quantity vector Q, and hence quadSub, uses
    // SeisSol's order (0=xx, 1=yy, 2=zz, 3=xy, 4=yz, 5=xz). The shear block is
    // therefore permuted, and the compliance has to be assembled in the *quantity*
    // order so that it can be contracted with the stress moments directly.
    //
    //   SeisSol quantity  0   1   2   3    4    5
    //   standard Voigt    1   2   3   6    4    5
    constexpr std::array<std::size_t, 6> ToVoigt{0, 1, 2, 5, 3, 4};

    const std::array<std::array<double, 6>, 6> voigtEntries{{
        {material.c11, material.c12, material.c13, material.c14, material.c15, material.c16},
        {material.c12, material.c22, material.c23, material.c24, material.c25, material.c26},
        {material.c13, material.c23, material.c33, material.c34, material.c35, material.c36},
        {material.c14, material.c24, material.c34, material.c44, material.c45, material.c46},
        {material.c15, material.c25, material.c35, material.c45, material.c55, material.c56},
        {material.c16, material.c26, material.c36, material.c46, material.c56, material.c66},
    }};

    Eigen::Matrix<double, 6, 6> stiffness{};
    for (std::size_t i = 0; i < 6; ++i) {
      for (std::size_t j = 0; j < 6; ++j) {
        stiffness(i, j) = voigtEntries[ToVoigt[i]][ToVoigt[j]];
      }
    }

    // A physically admissible stiffness is symmetric positive definite. A
    // degenerate one (e.g. a fluid, where the shear block vanishes) would make
    // the inverse blow up and surface much later as "Detected Inf/NaN in
    // energies", which is a confusing place to find out.
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(stiffness);
    const auto smallest = solver.eigenvalues().minCoeff();
    const auto largest = solver.eigenvalues().maxCoeff();
    if (!(smallest > 0) || smallest < 1e-12 * largest) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true)) {
        logWarning() << "Anisotropic stiffness tensor is singular or not positive definite"
                     << "(smallest eigenvalue" << smallest << ", largest" << largest
                     << "). The elastic energy of such cells is not well defined and will be"
                     << "reported as zero. Further occurrences are not reported.";
      }
      return AnisotropicEnergyData{};
    }

    const auto inverse = stiffness.inverse().eval();

    AnisotropicEnergyData data{};
    std::copy_n(inverse.data(), data.matS.size(), data.matS.begin());
    return data;
  }

  template <typename LinearViewT, typename QuadraticViewT>
  static std::array<double, EnergyCount>
      computeEnergies(const AnisotropicMaterial& material,
                      const AnisotropicMaterial::EnergyData& data,
                      const LinearViewT& linSub,
                      const QuadraticViewT& quadSub,
                      const Moments& /*moments*/,
                      std::size_t /*sim*/) {
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

    output[MomentumXIdx] = curMomentumX;
    output[MomentumYIdx] = curMomentumY;
    output[MomentumZIdx] = curMomentumZ;

    double curElasticEnergy = 0;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        curElasticEnergy += data.matS[i + 6 * j] * quadSub(i, j);
      }
    }

    output[ElasticEnergyIdx] = 0.5 * curElasticEnergy;
    output[ElasticKineticIdx] = curKineticEnergy;

    return output;
  }
};

} // namespace seissol::model
#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_ENERGY_H_
