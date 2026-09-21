// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ENERGYBASE_H_
#define SEISSOL_SRC_EQUATIONS_ENERGYBASE_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <type_traits>

namespace seissol::model {

/// Physical dimension of a reported energy, resolved to an SIUnit by the writer.
enum class EnergyUnit : std::uint8_t {
  Energy,   ///< J
  Power,    ///< W -- instantaneous rates, e.g. viscous dissipation
  Moment,   ///< Nm
  Momentum, ///< Ns
  Scalar,   ///< dimensionless
};

/**
 * Describes one reported energy.
 *
 * The descriptor list is the single source of truth: EnergyCount is derived from
 * it and the output indices are looked up by name, so the list and the values
 * computed for it cannot drift apart.
 *
 * `group` controls the terminal output; the CSV always contains every entry under
 * its `name`. Energies sharing a group are summed and reported on one line as
 *
 *     <groupLabel>: <total>, <shortLabel> x%, <shortLabel> y%, ...
 *
 * A group with a single member prints just the total. An energy without a group
 * is accumulated and written, but not printed on its own -- the writer may still
 * pick it up in a bespoke block, as it does for the momentum triple.
 */
struct EnergyDescriptor {
  std::string_view name;
  EnergyUnit unit = EnergyUnit::Energy;
  /// Energies sharing a group are reported together.
  std::string_view group;
  /// The heading for the group. Set on exactly one of its members.
  std::string_view groupLabel;
  /// This member's name in the percentage list. Omitted for single-member groups.
  std::string_view shortLabel;
};

namespace detail {

template <std::size_t A, std::size_t B>
constexpr std::array<EnergyDescriptor, A + B>
    concat(const std::array<EnergyDescriptor, A>& first,
           const std::array<EnergyDescriptor, B>& second) {
  std::array<EnergyDescriptor, A + B> result{};
  for (std::size_t i = 0; i < A; ++i) {
    result[i] = first[i];
  }
  for (std::size_t i = 0; i < B; ++i) {
    result[A + i] = second[i];
  }
  return result;
}

template <std::size_t A, std::size_t B, typename... Rest>
constexpr auto concat(const std::array<EnergyDescriptor, A>& first,
                      const std::array<EnergyDescriptor, B>& second,
                      const Rest&... rest) {
  return concat(concat(first, second), rest...);
}

/// Position of `name` in `list`, or list.size() if absent.
template <std::size_t N>
constexpr std::size_t indexOf(const std::array<EnergyDescriptor, N>& list, std::string_view name) {
  for (std::size_t i = 0; i < N; ++i) {
    if (list[i].name == name) {
      return i;
    }
  }
  return N;
}

/**
 * Compile-time sanity check for a descriptor list. Enforces that
 *  - every entry is named and names are unique (a truncated initializer list
 *    would leave empty names, which used to go unnoticed because EnergyCount was
 *    written out by hand),
 *  - every group carries exactly one groupLabel,
 *  - all members of a group share a unit, since they are summed,
 *  - a shortLabel is only given to a group member.
 */
template <std::size_t N>
constexpr bool descriptorsWellFormed(const std::array<EnergyDescriptor, N>& list) {
  for (std::size_t i = 0; i < N; ++i) {
    if (list[i].name.empty()) {
      return false;
    }
    for (std::size_t j = i + 1; j < N; ++j) {
      if (list[i].name == list[j].name) {
        return false;
      }
    }
    if (list[i].group.empty() && !(list[i].groupLabel.empty() && list[i].shortLabel.empty())) {
      return false;
    }
    if (!list[i].group.empty()) {
      std::size_t labels = 0;
      for (std::size_t j = 0; j < N; ++j) {
        if (list[j].group == list[i].group) {
          if (!list[j].groupLabel.empty()) {
            ++labels;
          }
          if (list[j].unit != list[i].unit) {
            return false;
          }
        }
      }
      if (labels != 1) {
        return false;
      }
    }
  }
  return true;
}

} // namespace detail

// Reusable fragments. Most equation systems report the same core set, so they
// compose their descriptor list from these rather than repeating it.

inline constexpr std::array MomentumEnergies{
    EnergyDescriptor{"momentumX", EnergyUnit::Momentum, {}, {}, {}},
    EnergyDescriptor{"momentumY", EnergyUnit::Momentum, {}, {}, {}},
    EnergyDescriptor{"momentumZ", EnergyUnit::Momentum, {}, {}, {}},
};

inline constexpr std::array AcousticEnergies{
    EnergyDescriptor{
        "acoustic_kinetic_energy", EnergyUnit::Energy, "acoustic", "Acoustic energy:", "kinematic"},
    EnergyDescriptor{"acoustic_energy", EnergyUnit::Energy, "acoustic", {}, "potential"},
};

inline constexpr std::array ElasticEnergies{
    EnergyDescriptor{
        "elastic_kinetic_energy", EnergyUnit::Energy, "elastic", "Elastic energy:", "kinematic"},
    EnergyDescriptor{"elastic_energy", EnergyUnit::Energy, "elastic", {}, "potential"},
};

/**
 * What a damaged cell reports.
 *
 * The free energy of the rheology is homogeneous of degree two in the strain,
 * so it is half the stress contracted with the strain -- and both are columns
 * of the transported tensor, which makes it a quadratic moment of that tensor.
 * Its difference from what the same strain would store in the undamaged solid
 * is what the damage has released. Of that difference the part linear in the
 * onset ratio is the damage times a quadratic form of the strain, which the
 * trilinear mass form integrates exactly; the part with the square root of the
 * second invariant is no polynomial and is reported as what remains. The two
 * need not share a sign: in compression the first invariant is negative, and
 * the root part with it.
 *
 * Kinetic and free energy are the two parts of the mechanical energy, and are
 * grouped as such; the released energy and its parts are not parts of it.
 */
inline constexpr std::array DamageEnergies{
    EnergyDescriptor{
        "damage_kinetic_energy", EnergyUnit::Energy, "damage", "Mechanical energy:", "kinematic"},
    EnergyDescriptor{"damage_free_energy", EnergyUnit::Energy, "damage", {}, "free"},
    EnergyDescriptor{"damage_undamaged_energy", EnergyUnit::Energy, {}, {}, {}},
    EnergyDescriptor{"damage_released_energy", EnergyUnit::Energy, {}, {}, {}},
    EnergyDescriptor{"damage_released_onset_energy", EnergyUnit::Energy, {}, {}, {}},
    EnergyDescriptor{"damage_released_root_energy", EnergyUnit::Energy, {}, {}, {}},
    EnergyDescriptor{"mean_damage", EnergyUnit::Scalar, {}, {}, {}},
    EnergyDescriptor{"mean_breakage", EnergyUnit::Scalar, {}, {}, {}},
};

/**
 * The Maxwell branch springs join the elastic group: they are part of the stored
 * potential energy, so reporting the kinetic/potential split without them would
 * understate the potential share.
 */
inline constexpr std::array ViscoelasticEnergies{
    EnergyDescriptor{"viscoelastic_energy", EnergyUnit::Energy, "elastic", {}, "viscoelastic"},
};

inline constexpr std::array ViscoacousticEnergies{
    EnergyDescriptor{"viscoacoustic_energy", EnergyUnit::Energy, "acoustic", {}, "viscoacoustic"},
};

inline constexpr std::array ViscoEnergies{
    EnergyDescriptor{"viscous_dissipation_rate",
                     EnergyUnit::Power,
                     "viscous_dissipation",
                     "Viscous dissipation rate:",
                     {}},
};

inline constexpr std::array DarcyEnergies{
    EnergyDescriptor{"darcy_dissipation_rate",
                     EnergyUnit::Power,
                     "darcy_dissipation",
                     "Darcy dissipation rate:",
                     {}},
};

/**
 * Specializations provide:
 *   Energies       -- constexpr std::array<EnergyDescriptor, N>, usually composed
 *                     from the fragments above via detail::concat
 *   EnergyCount    -- Energies.size(), never written out by hand
 *   EnergyData     -- per-cell data, see below
 *   Moments        -- per-cell moments beyond the ones EnergyOutput supplies
 *   initEnergyData -- builds EnergyData once, at setup
 *   computeMoments -- builds Moments for one cell
 *   computeEnergies-- evaluates the energies for one cell and one simulation
 *   WeightColumn   -- optional: the column of the transported tensor that weighs
 *                     the trilinear moment; see EnergyWeightColumn
 *
 * computeEnergies receives three moments of the transported tensor I over the
 * reference element: \\int I_J, \\int I_I I_J and \\int w I_I I_J, the last
 * with w the column a material names as WeightColumn. For a solver whose flux
 * is linear in its state, I is that state.
 *
 * Output positions are looked up by name with detail::indexOf rather than
 * hard-coded, so reordering the list cannot silently misplace a value.
 *
 * EnergyData is *precomputed, time-independent* per-cell data -- the anisotropic
 * compliance matrix is the motivating case. It is deliberately **not** an
 * accumulator: computeEnergies takes it by const reference and is a pure function
 * of the current state.
 *
 * Consequence for dissipative materials: report the instantaneous dissipation
 * *rate* (EnergyUnit::Power), not the time-integrated dissipated energy. The
 * integral is a state that would have to be advanced every time step, inside
 * kernels that run on the device, which is out of scope here; integrating the
 * reported rate in postprocessing gives the same quantity to the accuracy of the
 * output interval.
 */
template <typename MaterialT>
struct EnergyCompute;

/// A material that names no column has no use for the trilinear moment, and
/// the energy output does not form it.
inline constexpr std::size_t NoEnergyWeight = static_cast<std::size_t>(-1);

template <typename ComputeT, typename = void>
struct EnergyWeightColumn {
  static constexpr std::size_t Value = NoEnergyWeight;
};

template <typename ComputeT>
struct EnergyWeightColumn<ComputeT, std::void_t<decltype(ComputeT::WeightColumn)>> {
  static constexpr std::size_t Value = ComputeT::WeightColumn;
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ENERGYBASE_H_
