// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MODEL_CONCEPTS_H_
#define SEISSOL_SRC_MODEL_CONCEPTS_H_

#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/Solver.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"
#include "Model/CommonDatastructures.h"

#include <concepts>

namespace seissol {

template <typename T>
concept SolverC =
    requires {
      typename T::TimeKernelT;
      typename T::SpacetimeKernelT;
      typename T::LocalKernelT;
      typename T::NeighborKernelT;

      { T::DerivativesSize } -> std::same_as<const std::size_t&>;
      { T::BuffersSize } -> std::same_as<const std::size_t&>;
    } && std::derived_from<typename T::TimeKernelT, kernels::TimeKernel> &&
    std::derived_from<typename T::SpacetimeKernelT, kernels::SpacetimeKernel> &&
    std::derived_from<typename T::LocalKernelT, kernels::LocalKernel> &&
    std::derived_from<typename T::NeighborKernelT, kernels::NeighborKernel>;

template <typename T>
concept ComponentC = requires {
  { T::NumQuantities } -> std::same_as<const std::size_t&>;
  { T::Text } -> std::same_as<const std::string&>;
  { T::Quantities } -> std::same_as<const std::array<std::string, T::NumQuantities>&>;
  { T::Parameters } -> std::same_as<const std::size_t&>;
};

template <typename T>
concept MaterialC = ComponentC<T> && requires {
  { T::NumberPerMechanism } -> std::same_as<const std::size_t&>;
  { T::TractionQuantities } -> std::same_as<const std::size_t&>;
  { T::Mechanisms } -> std::same_as<const std::size_t&>;
  { T::Type } -> std::same_as<const model::MaterialType&>;

  { T::SupportsDR } -> std::same_as<const bool&>;
  { T::SupportsLTS } -> std::same_as<const bool&>;

  typename T::LocalSpecificData;
  typename T::NeighborSpecificData;
  typename T::Solver;
} && SolverC<typename T::Solver>;

template <typename T>
concept ConfigC = requires {
  { T::ConvergenceOrder } -> std::same_as<const std::size_t&>;
  { T::RelaxationMechanisms } -> std::same_as<const std::size_t&>;
  { T::MaterialType } -> std::same_as<const model::MaterialType&>;
  { T::Precision } -> std::same_as<const RealType&>;
  { T::ViscoMode } -> std::same_as<const ViscoImplementation&>;
  { T::DRQuadRule } -> std::same_as<const DRQuadRuleType&>;
  { T::NumSimulations } -> std::same_as<const std::size_t&>;
};

} // namespace seissol
#endif // SEISSOL_SRC_MODEL_CONCEPTS_H_
