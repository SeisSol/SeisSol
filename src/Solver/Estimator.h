// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: David Schneller

#ifndef SEISSOL_SRC_SOLVER_ESTIMATOR_H_
#define SEISSOL_SRC_SOLVER_ESTIMATOR_H_

#include "Common/ConfigRegistry.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace seissol::solver {

/// The time the local kernel of `config` takes on this rank for a fixed amount of work.
auto miniSeisSol(ConfigId config) -> double;

/// How the work of a cell depends on its configuration: the cost of a cell of each configuration in
/// `configs` relative to a cell of `reference`, indexed by the id of the configuration; 1 for the
/// configurations built that are not in `configs`.
///
/// With `measure`, every rank times the time, local and neighbor kernels of each configuration with
/// the proxy, and all ranks take the median of the ratios over the ranks; collective then.
/// Otherwise, the ratios are those of the hardware FLOPs the kernels count, weighted with the size
/// of a real; that sees nothing of the machine.
auto configCostFactors(const std::vector<ConfigId>& configs, ConfigId reference, bool measure)
    -> std::vector<double>;

auto hostDeviceSwitch() -> int;

/**
 * What a step of the time, local and neighbor kernels of a layer costs over the number of its
 * cells: measured at a few sizes, linear in between and, below the smallest and beyond the largest
 * size, along the first and the last stretch (not below zero). A layer of no cells costs nothing.
 * The measured times become non-decreasing, so that noise cannot make a larger layer cheaper than a
 * smaller one.
 */
class LayerCost {
  public:
  LayerCost() = default;

  /// `points`: numbers of cells, increasing, with the time of a step of a layer of that many
  explicit LayerCost(std::vector<std::pair<std::size_t, double>> points);

  /// the time of a step of a layer of `cells` cells
  [[nodiscard]] double operator()(std::size_t cells) const;

  /**
   * What one layer of `first + second` cells saves over a layer of `first` and one of `second`
   * cells, per step: the fixed part of the second layer, and what the smaller layers do not fill of
   * the device or the threads; nothing if one of them is empty.
   */
  [[nodiscard]] double combinedSaving(std::size_t first, std::size_t second) const;

  [[nodiscard]] const std::vector<std::pair<std::size_t, double>>& points() const {
    return points_;
  }

  private:
  std::vector<std::pair<std::size_t, double>> points_;
};

/**
 * Times a step of the time, local and neighbor kernels of a layer of `config` with the proxy, with
 * 8 cells and twice as many each time up to 16384, or until a step takes longer than 20 ms on all
 * ranks; on every rank (on the device in device builds, with the threads of the rank otherwise).
 * Each size is the median over the ranks: the same on all of them. Collective; takes a few
 * seconds (some 20 ms per size and run, the setup of each run, a warm-up of 0.3 s).
 */
auto measureLayerCost(ConfigId config) -> LayerCost;

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_ESTIMATOR_H_
