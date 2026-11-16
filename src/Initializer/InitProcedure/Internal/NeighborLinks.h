// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_NEIGHBORLINKS_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_NEIGHBORLINKS_H_

#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Tree/Colormap.h"

#include <cstddef>
#include <limits>
#include <map>
#include <vector>

namespace seissol::initializer::internal {

/// How long a message from a neighbor rank takes: a latency, and a bandwidth beyond it.
struct NeighborLink {
  double latency{0};
  double bytesPerSecond{std::numeric_limits<double>::infinity()};

  /// the time a message of `bytes` takes
  [[nodiscard]] double time(std::size_t bytes) const {
    return latency + static_cast<double>(bytes) / bytesPerSecond;
  }
};

/**
 * Measures the links from the ranks that the halo of `layout` exchanges cells with, as the halo
 * exchange uses them: with messages to and from all of them at once, once small and once large
 * (1 MiB), the best of three runs each; in host memory. Collective.
 */
std::map<int, NeighborLink> measureNeighborLinks(const MeshLayout& layout);

/**
 * Per color, for the copy layers: how long the halo cells take to arrive that the cluster of the
 * layer waits for, as `links` estimate it: the slowest of the ranks they come from, with
 * `bytesPerCell` (per color) for each cell of a ghost layer. 0 for the other colors and for copy
 * layers that wait for no cells.
 */
std::vector<double> haloTimes(const MeshLayout& layout,
                              const LTSColorMap& colormap,
                              const std::map<int, NeighborLink>& links,
                              const std::vector<std::size_t>& bytesPerCell);

} // namespace seissol::initializer::internal

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_NEIGHBORLINKS_H_
