// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MODEL_CONCEPTS_H_
#define SEISSOL_SRC_MODEL_CONCEPTS_H_

#include "Common/ConfigValue.h"
#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"
#include "Model/CommonDatastructures.h"
#include "Model/MaterialType.h"
#include "Model/Quantities.h"

#include <array>
#include <concepts>
#include <cstddef>
#include <ranges>
#include <string>

namespace seissol {

namespace internal {
/// An array of quantity groups, as a material lays out its quantities in them.
template <typename T>
concept QuantityGroupsC = std::ranges::random_access_range<T> &&
                          std::same_as<std::ranges::range_value_t<T>, model::QuantityGroup>;
} // namespace internal

/// A configuration, as `ConfigOf` defines them: the constants the code that is generic over
/// configurations reads.
template <typename T>
concept ConfigC = requires {
  { T::Value } -> std::same_as<const ConfigValue&>;
  { T::ConvergenceOrder } -> std::same_as<const std::size_t&>;
  { T::RelaxationMechanisms } -> std::same_as<const std::size_t&>;
  { T::MaterialType } -> std::same_as<const model::MaterialType&>;
  { T::Precision } -> std::same_as<const RealType&>;
  { T::Solver } -> std::same_as<const SolverType&>;
  { T::DRQuadRule } -> std::same_as<const DRQuadRuleType&>;
  { T::NumSimulations } -> std::same_as<const std::size_t&>;
};

/// A material model: a `model::Material` with the static properties that the code that is generic
/// over materials reads. Its quantity groups for the rotations to and from a face depend on the
/// solver; Model/Common.h checks them against the generated layout.
template <typename T>
concept MaterialC =
    std::derived_from<T, model::Material> &&
    requires {
      { T::Type } -> std::same_as<const model::MaterialType&>;
      { T::Text } -> std::same_as<const std::string&>;
      { T::NumQuantities } -> std::same_as<const std::size_t&>;
      { T::NumElasticQuantities } -> std::same_as<const std::size_t&>;
      { T::NumberPerMechanism } -> std::same_as<const std::size_t&>;
      { T::Mechanisms } -> std::same_as<const std::size_t&>;
      { T::Parameters } -> std::same_as<const std::size_t&>;
      { T::VelocityOffset } -> std::same_as<const std::size_t&>;
      { T::TractionComponents } -> std::same_as<const std::size_t&>;
      // the names of the quantities, without the memory variables of the mechanisms
      { T::Quantities } -> std::same_as<const std::array<std::string, T::NumElasticQuantities>&>;
      { T::SupportsDR } -> std::same_as<const bool&>;
      { T::SupportsPlasticity } -> std::same_as<const bool&>;
      { T::SupportsLTS } -> std::same_as<const bool&>;
      { T::SupportsEnergy } -> std::same_as<const bool&>;
      typename T::LocalSpecificData;
      typename T::NeighborSpecificData;
    } && std::derived_from<typename T::RiemannMaterial, model::Material> &&
    internal::QuantityGroupsC<decltype(T::PrimaryGroups)> &&
    std::same_as<std::ranges::range_value_t<decltype(T::StiffSourceRows)>, model::StiffSourceRow>;

/// A solver, as `kernels::SolverOf` selects one for the configuration `Cfg`: its kernels and the
/// properties of its data.
template <typename T, typename Cfg>
concept SolverC =
    requires {
      typename T::LocalData;
      typename T::NeighborData;
      typename T::template TimeBasis<Real<Cfg>>;
      { T::FusedMechanisms } -> std::same_as<const bool&>;
      { T::IntegralsSize } -> std::same_as<const std::size_t&>;
      { T::DerivativesSize } -> std::same_as<const std::size_t&>;
      { T::implementsFaceType(FaceType::Regular) } -> std::same_as<FaceTypeSupport>;
    } && std::derived_from<typename T::TimeKernelT, kernels::TimeKernel<Cfg>> &&
    std::derived_from<typename T::SpacetimeKernelT, kernels::SpacetimeKernel<Cfg>> &&
    std::derived_from<typename T::LocalKernelT, kernels::LocalKernel<Cfg>> &&
    std::derived_from<typename T::NeighborKernelT, kernels::NeighborKernel<Cfg>>;

} // namespace seissol
#endif // SEISSOL_SRC_MODEL_CONCEPTS_H_
