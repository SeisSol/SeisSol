// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/TimeStepping/TimestepHistogram.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <mpi.h>
#include <utility>
#include <vector>

namespace seissol::initializer {

TimestepHistogram::TimestepHistogram(std::vector<double> cumulative)
    : cumulative_(std::move(cumulative)) {}

TimestepHistogram TimestepHistogram::fromCells(const std::vector<double>& cellTimesteps,
                                               const std::vector<int>& cellCosts,
                                               double baseTimestep,
                                               std::size_t maxIndex) {
  assert(cellTimesteps.size() == cellCosts.size());
  assert(baseTimestep > 0.0);

  std::vector<double> bins(maxIndex + 1, 0.0);
  const auto maxIndexAsDouble = static_cast<double>(maxIndex);

  for (std::size_t cell = 0; cell < cellTimesteps.size(); ++cell) {
    const auto scaled = cellTimesteps[cell] / baseTimestep;
    // clamp before the cast: a timestep far above the ladder would otherwise wrap
    const auto index =
        scaled >= maxIndexAsDouble ? maxIndex : static_cast<std::size_t>(std::max(scaled, 0.0));
    bins[index] += cellCosts[cell];
  }

  std::vector<double> cumulative(maxIndex + 2, 0.0);
  for (std::size_t index = 0; index <= maxIndex; ++index) {
    cumulative[index + 1] = cumulative[index] + bins[index];
  }
  return TimestepHistogram(std::move(cumulative));
}

void TimestepHistogram::reduce(MPI_Comm comm) {
  if (cumulative_.empty()) {
    return;
  }
  MPI_Allreduce(MPI_IN_PLACE,
                cumulative_.data(),
                static_cast<int>(cumulative_.size()),
                MPI_DOUBLE,
                MPI_SUM,
                comm);
}

double TimestepHistogram::weightBelow(std::size_t index) const {
  if (cumulative_.empty()) {
    return 0.0;
  }
  return cumulative_[std::min(index, cumulative_.size() - 1)];
}

} // namespace seissol::initializer
