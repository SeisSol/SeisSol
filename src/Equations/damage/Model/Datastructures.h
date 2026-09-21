// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DATASTRUCTURES_H_

#include "Kernels/NonLinearCK/Data.h"
#include "Kernels/SolverSelector.h"
#include "Model/CommonDatastructures.h"
#include "Model/Quantities.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace seissol::model {

/// Continuum damage-breakage rheology.
///
/// The state is a strain rather than a stress, which is what lets the two
/// internal variables enter the constitutive law instead of the flux: alpha
/// degrades the moduli, B blends the solid branch into a granular one. Both
/// evolve through source terms local to the cell, so neither takes part in the
/// face coupling.
///
/// Only the quantities that vary per cell live here. The coefficients of the
/// granular branch and the two rates are properties of the model rather than
/// of a point in the mesh, and come from the parameter file.
struct DamageMaterial : public Material {
  static constexpr std::size_t NumQuantities = 11;
  static constexpr std::size_t NumElasticQuantities = 9;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Mechanisms = 0;
  static constexpr MaterialType Type = MaterialType::Damage;
  static inline const std::string Text = "damage";
  static inline const std::array<std::string, NumQuantities> Quantities = {"eps_xx",
                                                                           "eps_yy",
                                                                           "eps_zz",
                                                                           "eps_xy",
                                                                           "eps_yz",
                                                                           "eps_xz",
                                                                           "v1",
                                                                           "v2",
                                                                           "v3",
                                                                           "alpha",
                                                                           "breakage"};

  using Solver = kernels::SolverSelector<Config::Solver>::Type;

  static constexpr auto PrimaryGroups = DamageQuantities;
  static constexpr auto TransportGroups = DamageTransportQuantities;
  static constexpr auto RotationGroups = PrimaryGroups;
  static constexpr auto InverseRotationGroups = PrimaryGroups;

  static constexpr std::size_t VelocityOffset = roleOffset(PrimaryGroups, FaceRole::Velocity);
  static constexpr std::size_t TractionComponents = roleExtent(PrimaryGroups, FaceRole::Traction);

  /// rho, plus the five moduli and rates and the six components of the
  /// initial strain.
  static constexpr std::size_t Parameters = 11 + Material::Parameters;

  /// A fault reads the traction out of what a cell transports rather than out
  /// of its state, and the impedance it scales with is formed at the node from
  /// the state that is there.
  static constexpr bool SupportsDR = true;
  // What a cell transports beyond its state is projected into an expansion of
  // its own, so a neighbour on a coarser cluster reconstructs a subinterval of
  // the stress from that -- rather than rebuilding the stress out of this
  // cell's material, which is a question this cell answers and not its
  // neighbour. Whether the two agree to the scheme's order is what a run
  // against global time stepping says, and that has not been run.
  static constexpr bool SupportsLTS = true;
  static constexpr bool SupportsEnergy = true;

  using LocalSpecificData = kernels::solver::nonlinearck::NonLinearLocalData;
  using NeighborSpecificData = kernels::solver::nonlinearck::NonLinearNeighborData;

  using EnergyData = std::monostate;

  /// Undamaged Lame parameters. The effective ones are not stored: they depend
  /// on the strain state through the invariant ratio, so they are a property
  /// of a point and not of a cell, and the kernels form them where they are
  /// used.
  double lambda0{};
  double mu0{};
  /// Modulus coupling damage to the strain invariants.
  double gammaR{};
  /// Onset of damage growth, as a value of the invariant ratio.
  double xi0{};
  /// Rate at which damage accumulates above the onset.
  double damageRate{};
  /// Strain the cell already carries before the simulation starts, in Voigt
  /// order. It enters every invariant, so it is not an initial condition on Q
  /// but a property of the material -- and one that varies with depth in every
  /// scenario that loads a fault, which is why the six components are queried
  /// one by one rather than held as an array: the map binds a member per
  /// supplied parameter.
  /// Properties of the model rather than of a point in the mesh, so they come
  /// from the parameter file rather than from easi, and they arrive through
  /// initialize() like every other parameter that is read there.
  double breakageRate{};
  double healingRate{};
  /// Width of the smoothed step from damage to breakage. It divides, so it is
  /// not allowed to be zero even where breakage is switched off.
  double betaAlpha{1.0};
  std::array<double, 4> aB{};

  double epsInitXX{};
  double epsInitYY{};
  double epsInitZZ{};
  double epsInitXY{};
  double epsInitYZ{};
  double epsInitXZ{};

  static const std::unordered_map<std::string, double DamageMaterial::*> ParameterMap;

  [[nodiscard]] double getLambdaBar() const override { return lambda0; }

  [[nodiscard]] double getMuBar() const override { return mu0; }

  /// Speeds of the undamaged material, which is what the moduli describe.
  [[nodiscard]] double getPWaveSpeed() const override {
    return std::sqrt((lambda0 + 2.0 * mu0) / rho);
  }

  [[nodiscard]] double getSWaveSpeed() const override { return std::sqrt(mu0 / rho); }

  /// Fastest wave the material can carry, over every state it can be in.
  ///
  /// What a wave travels at follows from the tangent of the stress, not from
  /// its ratio to the strain. With n the strain normalised in the Frobenius
  /// norm, xi = tr(n) and g = gammaR alpha, differentiating the solid branch
  /// gives
  ///
  ///   C = lambda0 d(x)d + (2 mu0 - 2 g xi0 - g xi) Isym
  ///         - g (d(x)n + n(x)d) + g xi n(x)n,
  ///
  /// whose last two groups no pair of Lame parameters can express. In a
  /// direction e the acoustic tensor is m Id + S, with m = mu0 - g xi0 -
  /// g xi / 2 and S supported on the plane spanned by e and n.e. Its largest
  /// eigenvalue is maximal where n.e is parallel to e, and there, with
  /// a = e.n.e,
  ///
  ///   rho c^2 = lambda0 + 2 mu0 + g (xi (a^2 - 1) - 2 a - 2 xi0).
  ///
  /// Since a^2 <= 1 this grows as xi falls, and |n| = 1 holds xi above
  /// a - sqrt(2 (1 - a^2)); that floor reaches the range limit -sqrt(3) at
  /// a = -1 / sqrt(3), and only there. The excess over the undamaged
  /// stiffness is therefore 4 / sqrt(3) - 2 xi0 per unit of gammaR, all of it
  /// carried at alpha = 1, and negative xi0 is the case where damage only
  /// softens and the undamaged speed already bounds everything.
  ///
  /// The size of the strain drops out: the tangent depends on its direction
  /// alone, so this is a constant of the material rather than something a
  /// cell has to be asked for. It bounds the solid branch; a granular branch
  /// with non-zero aB would add a term of its own.
  /// The stiffness of the undamaged solid, which is what an interface between
  /// two cells is judged by at setup. The damaged one is a function of the
  /// state and has no place in a tensor that is filled once.
  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {
    auto view = seissol::general::init::stiffnessTensor::view::create(fullTensor.data());
    view.setZero();
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        view(i, i, j, j) = lambda0;
      }
      view(i, i, i, i) = lambda0 + 2.0 * mu0;
      for (std::size_t j = 0; j < 3; ++j) {
        if (i != j) {
          view(i, j, i, j) = mu0;
          view(i, j, j, i) = mu0;
        }
      }
    }
  }

  /// Stiffest the granular branch can be, over every strain direction.
  ///
  /// Its tangent has the same four groups the solid one has, and in both of
  /// them the directional group is minus the invariant ratio times the mixed
  /// one. So the acoustic tensor in a direction e comes out as
  ///
  ///   rho c^2 = L + M + X (2 a - xi a^2),   a = e.n.e,
  ///
  /// with L the weight on d(x)d, M the one on Isym and X the one on the mixed
  /// group -- for the granular branch L = P'', M = 2P - xi P' and X = P' - xi
  /// P''. The dependence on a is a quadratic and is maximised in closed form;
  /// what is left is the invariant ratio, which the coefficients are a cubic in
  /// at worst, and which is scanned. The constraint that ties a to xi is not
  /// imposed, and neither is a negative bound kept: both err upwards, which
  /// costs a timestep rather than a run.
  ///
  /// With aB1 = aB3 = 0 the mixed group is zero and this is 2 aB0 + 2 aB2, the
  /// isotropic pair of a granular branch that has one.
  [[nodiscard]] double granularStiffness() const {
    constexpr double RatioLimit = 1.7320508075688772; // sqrt(3)
    constexpr std::size_t Samples = 257;
    double stiffest = 0.0;
    for (std::size_t sample = 0; sample < Samples; ++sample) {
      const double ratio =
          -RatioLimit + 2.0 * RatioLimit * static_cast<double>(sample) / (Samples - 1);
      const double volumetric = 2.0 * aB[2] + 6.0 * aB[3] * ratio;
      const double isotropic = 2.0 * aB[0] + aB[1] * ratio - aB[3] * ratio * ratio * ratio;
      const double mixed = aB[1] - 3.0 * aB[3] * ratio * ratio;
      double directional = std::max(mixed * (2.0 - ratio), -mixed * (2.0 + ratio));
      if (mixed * ratio > 0.0 && std::abs(ratio) >= 1.0) {
        directional = std::max(directional, mixed / ratio);
      }
      stiffest = std::max(stiffest, volumetric + isotropic + directional);
    }
    return stiffest;
  }

  [[nodiscard]] double getMaxWaveSpeed() const override {
    constexpr double TangentExcess = 2.3094010767585034; // 4 / sqrt(3)
    const double damage = std::max(gammaR * (TangentExcess - 2.0 * xi0), 0.0);
    const double solid = lambda0 + 2.0 * mu0 + damage;
    // The stress is a convex blend of the two branches, so the acoustic tensor
    // is, and so the larger of the two bounds is one for the blend. Taking the
    // solid branch alone was a bound only for a medium that never breaks.
    return std::sqrt(std::max(solid, granularStiffness()) / rho);
  }

  DamageMaterial() = default;
  explicit DamageMaterial(const std::vector<double>& materialValues)
      : Material(materialValues), lambda0(materialValues.at(1)), mu0(materialValues.at(2)),
        gammaR(materialValues.at(3)), xi0(materialValues.at(4)), damageRate(materialValues.at(5)),
        epsInitXX(materialValues.at(6)), epsInitYY(materialValues.at(7)),
        epsInitZZ(materialValues.at(8)), epsInitXY(materialValues.at(9)),
        epsInitYZ(materialValues.at(10)), epsInitXZ(materialValues.at(11)) {}

  ~DamageMaterial() override = default;

  void initialize(const initializer::parameters::ModelParameters& parameters) override {
    const auto& damage = parameters.damageParameters;
    breakageRate = damage.breakageRate;
    healingRate = damage.healingRate;
    betaAlpha = damage.betaAlpha;
    aB = damage.granular;
  }

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};

inline const std::unordered_map<std::string, double DamageMaterial::*> DamageMaterial::ParameterMap{
    {"rho", &DamageMaterial::rho},
    {"lambda0", &DamageMaterial::lambda0},
    {"mu0", &DamageMaterial::mu0},
    {"gammaR", &DamageMaterial::gammaR},
    {"xi0", &DamageMaterial::xi0},
    {"Cd", &DamageMaterial::damageRate},
    {"eps_xx0", &DamageMaterial::epsInitXX},
    {"eps_yy0", &DamageMaterial::epsInitYY},
    {"eps_zz0", &DamageMaterial::epsInitZZ},
    {"eps_xy0", &DamageMaterial::epsInitXY},
    {"eps_yz0", &DamageMaterial::epsInitYZ},
    {"eps_xz0", &DamageMaterial::epsInitXZ},
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_DATASTRUCTURES_H_
