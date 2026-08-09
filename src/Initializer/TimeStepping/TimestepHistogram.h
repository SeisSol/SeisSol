// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_TIMESTEPHISTOGRAM_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_TIMESTEPHISTOGRAM_H_

#include <cstddef>
#include <mpi.h>
#include <vector>

namespace seissol::initializer {

/// Cell cost accumulated over `floor(cellTimestep / baseTimestep)`, stored cumulatively.
///
/// Unlike ClusterHistogram this is indexed by the update factor itself rather than by a
/// cluster id, so it does not presuppose a ladder. Any ladder built on the same base
/// timestep can be scored from it: cluster k holds exactly the cells with index in
/// `[updateFactor(k), updateFactor(k+1))`, which is one subtraction. That is what lets a
/// search enumerate ladders without touching the mesh or communicating again.
///
/// # Binning is not bit-identical to ClusterLadder
///
/// This compares `cellTimestep / baseTimestep` against integers, whereas ClusterLadder
/// accumulates its boundaries as doubles. The two can disagree for a cell sitting exactly on
/// a boundary. That is harmless as long as this type is only used to *rank* candidates and
/// the winner is then realized through ClusterLadder, which is how the search is meant to
/// work -- but it does mean a predicted cost can differ marginally from the realized one.
class TimestepHistogram {
  public:
  TimestepHistogram() = default;

  /// Bins cells, folding everything at or above `maxIndex` into the topmost bin.
  static TimestepHistogram fromCells(const std::vector<double>& cellTimesteps,
                                     const std::vector<int>& cellCosts,
                                     double baseTimestep,
                                     std::size_t maxIndex);

  /// Sums across all ranks in place. Cumulative sums add elementwise, so this is the same as
  /// reducing the underlying histogram.
  void reduce(MPI_Comm comm);

  /// Largest index that has its own bin; everything above is folded into it.
  [[nodiscard]] std::size_t maxIndex() const {
    return cumulative_.size() < 2 ? 0 : cumulative_.size() - 2;
  }

  /// Total weight of all cells whose index is strictly below `index`.
  [[nodiscard]] double weightBelow(std::size_t index) const;

  /// Total weight of cells with index in `[begin, end)`.
  [[nodiscard]] double weightIn(std::size_t begin, std::size_t end) const {
    return weightBelow(end) - weightBelow(begin);
  }

  [[nodiscard]] double totalWeight() const {
    return cumulative_.empty() ? 0.0 : cumulative_.back();
  }

  [[nodiscard]] const std::vector<double>& cumulative() const { return cumulative_; }

  private:
  explicit TimestepHistogram(std::vector<double> cumulative);

  /// cumulative_[i] is the weight of all cells with index < i; size is maxIndex + 2.
  std::vector<double> cumulative_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_TIMESTEPHISTOGRAM_H_
