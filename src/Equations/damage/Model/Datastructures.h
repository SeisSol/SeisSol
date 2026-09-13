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

#include <array>
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
  static constexpr auto RotationGroups = PrimaryGroups;
  static constexpr auto InverseRotationGroups = PrimaryGroups;

  static constexpr std::size_t VelocityOffset = roleOffset(PrimaryGroups, FaceRole::Velocity);
  static constexpr std::size_t TractionComponents = roleExtent(PrimaryGroups, FaceRole::Traction);

  /// rho, plus the five moduli and rates and the six components of the
  /// initial strain.
  static constexpr std::size_t Parameters = 11 + Material::Parameters;

  /// Dynamic rupture needs a mechanical traction, which is a derived quantity
  /// here and not a state variable. Until that is settled the material says so
  /// rather than producing something that looks like a traction and is not.
  static constexpr bool SupportsDR = false;
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
  double epsInitXX{};
  double epsInitYY{};
  double epsInitZZ{};
  double epsInitXY{};
  double epsInitYZ{};
  double epsInitXZ{};

  static const std::unordered_map<std::string, double DamageMaterial::*> ParameterMap;

  [[nodiscard]] double getLambdaBar() const override { return lambda0; }

  [[nodiscard]] double getMuBar() const override { return mu0; }

  DamageMaterial() = default;
  explicit DamageMaterial(const std::vector<double>& materialValues)
      : Material(materialValues), lambda0(materialValues.at(1)), mu0(materialValues.at(2)),
        gammaR(materialValues.at(3)), xi0(materialValues.at(4)), damageRate(materialValues.at(5)),
        epsInitXX(materialValues.at(6)), epsInitYY(materialValues.at(7)),
        epsInitZZ(materialValues.at(8)), epsInitXY(materialValues.at(9)),
        epsInitYZ(materialValues.at(10)), epsInitXZ(materialValues.at(11)) {}

  ~DamageMaterial() override = default;

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
