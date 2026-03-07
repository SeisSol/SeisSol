// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Initializer/TimeStepping/ClusterLayout.h"

#include "Geometry/MeshReader.h"
#include "Monitoring/Unit.h"
#include "Numerical/StableSum.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mpi.h>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

ClusterLayout ClusterLayout::fromMesh(const std::vector<std::uint64_t>& rates,
                                      const geometry::MeshReader& mesh,
                                      double wiggle,
                                      bool infoprint) {

  if (infoprint) {
    std::size_t cellCount = mesh.getElements().size();
    const auto summary = statistics::parallelSummary(cellCount);

    cellCount = seissol::Mpi::mpi.allreduce(cellCount, MPI_SUM);

    logInfo() << "Cell count:" << formatInteger(cellCount).c_str() << "(per rank:" << summary.mean
              << "±" << summary.std << "; range: [" << summary.min << ";" << summary.max << "])";
  }

  std::uint64_t maxLtsId = 0;
  double minimumTimestep = std::numeric_limits<double>::max();
  for (const auto& element : mesh.getElements()) {
    maxLtsId = std::max(maxLtsId, static_cast<std::uint64_t>(element.clusterId));
    minimumTimestep = std::min(minimumTimestep, element.timestep);
  }

  maxLtsId = seissol::Mpi::mpi.allreduce(maxLtsId, MPI_MAX);
  minimumTimestep = seissol::Mpi::mpi.allreduce(minimumTimestep, MPI_MIN);

  if (wiggle == 1) {
    if (infoprint) {
      logInfo() << "Minimum timestep:" << seissol::UnitTime.formatPrefix(minimumTimestep).c_str();
    }
  } else {
    if (infoprint) {
      logInfo() << "Minimum timestep (pre-wiggle):"
                << seissol::UnitTime.formatPrefix(minimumTimestep).c_str();
    }

    // apply wiggle here
    minimumTimestep *= wiggle;
    if (infoprint) {
      logInfo() << "Minimum timestep (with wiggle" << wiggle
                << "):" << seissol::UnitTime.formatPrefix(minimumTimestep).c_str();
    }
  }
  ClusterLayout layout(rates, minimumTimestep, maxLtsId + 1);

  if (infoprint) {
    std::vector<std::uint64_t> clusters(maxLtsId + 1);
    std::vector<std::uint64_t> clustersDR(maxLtsId + 1);
    const double timefactorGTS = minimumTimestep * mesh.getElements().size();
    numerical::StableAccumulator<double> timefactorELTS;
    numerical::StableAccumulator<double> timefactorCLTS;
    for (const auto& element : mesh.getElements()) {
      ++clusters[static_cast<std::size_t>(element.clusterId)];
      for (int i = 0; i < 4; ++i) {
        if (element.boundaries[i] == 3) {
          ++clustersDR[static_cast<std::size_t>(element.clusterId)];
        }
      }

      const auto elts = element.timestep * wiggle;
      const auto clts = layout.timestepRate(static_cast<std::size_t>(element.clusterId));

      timefactorELTS += elts;
      timefactorCLTS += clts;
    }

    seissol::Mpi::mpi.allreduceContainer(clusters, MPI_SUM);
    seissol::Mpi::mpi.allreduceContainer(clustersDR, MPI_SUM);

    std::array<double, 3> timefactors{
        timefactorGTS, timefactorELTS.result(), timefactorCLTS.result()};
    seissol::Mpi::mpi.allreduceContainer(timefactors, MPI_SUM);

    const auto eltsSpeedup = timefactors[1] / timefactors[0];
    const auto cltsSpeedup = timefactors[2] / timefactors[0];
    logInfo() << "Theoretical speedup to GTS:" << eltsSpeedup << "elementwise LTS;" << cltsSpeedup
              << "clustered LTS (current setup)";

    logInfo() << "Cluster histogram (cells, DR faces):";
    for (std::size_t i = 0; i <= maxLtsId; ++i) {
      // NOTE: we count the DR faces twice; hence divide by 2 here
      // (also we now ignore copy-layer DR faces in the histogram unlike @:1.3.2; see the actual
      // clustering histogram for that instead)
      logInfo() << i << ":" << clusters[i] << "," << (clustersDR[i] / 2);
    }
  }

  return layout;
}

} // namespace seissol::initializer
