// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERCOSTMODEL_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERCOSTMODEL_H_

#include <algorithm>
#include <cstdint>

namespace seissol::initializer {

/// What one time cluster costs, per unit of simulated time and up to a factor of the base
/// timestep.
///
/// SeisSol has always counted cell updates: a cluster carrying weight W that updates every
/// n-th base timestep costs W/n. That is exactly this model with `launchCost == 0` and
/// `fillThreshold == 0` -- not an approximation of it, the same expression, so a search
/// configured that way reproduces the old objective bit for bit.
///
/// The two extra terms cover what update counting misses:
///   launchCost     paid per cluster update no matter how little work it carries (kernel
///                  launches, the scheduler round trip, the communication handshake)
///   fillThreshold  work below which an update takes as long as one carrying exactly that
///                  much work, because the device never saturates
///
/// Both make a cluster that holds almost nothing expensive, which is what makes splitting a
/// ladder ever stop being free.
struct ClusterCostModel {
  double updateCost{1.0};
  double launchCost{0.0};
  double fillThreshold{0.0};

  /// Contribution of one cluster, still to be divided by the base timestep.
  [[nodiscard]] double clusterTerm(double weight, std::uint64_t updateFactor) const {
    return (launchCost + updateCost * std::max(weight, fillThreshold)) /
           static_cast<double>(updateFactor);
  }

  /// True if this is the pure update-count objective.
  [[nodiscard]] bool isUpdateCount() const {
    return launchCost == 0.0 && fillThreshold == 0.0 && updateCost == 1.0;
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERCOSTMODEL_H_
