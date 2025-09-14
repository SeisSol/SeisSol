// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include <Geometry/MeshReader.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Monitoring/Unit.h>
#include <Numerical/StableSum.h>
#include <Parallel/MPI.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utils/logger.h>
#include <vector>

#include <mpi.h>

namespace seissol::initializer {

ClusterLayout ClusterLayout::fromMesh(const std::vector<std::uint64_t>& rates,
                                      const geometry::MeshReader& mesh,
                                      double wiggle,
                                      bool infoprint) {
  std::uint64_t maxLtsId = 0;
  double minimumTimestep = std::numeric_limits<double>::max();
  for (const auto& element : mesh.getElements()) {
    maxLtsId = std::max(maxLtsId, static_cast<std::uint64_t>(element.clusterId));
    minimumTimestep = std::min(minimumTimestep, element.timestep);
  }
  MPI_Allreduce(
      MPI_IN_PLACE, &maxLtsId, 1, MPI::castToMpiType<std::uint64_t>(), MPI_MAX, MPI::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE, &minimumTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI::mpi.comm());
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
    MPI_Allreduce(MPI_IN_PLACE,
                  clusters.data(),
                  clusters.size(),
                  MPI::castToMpiType<std::uint64_t>(),
                  MPI_SUM,
                  MPI::mpi.comm());
    MPI_Allreduce(MPI_IN_PLACE,
                  clustersDR.data(),
                  clustersDR.size(),
                  MPI::castToMpiType<std::uint64_t>(),
                  MPI_SUM,
                  MPI::mpi.comm());
    std::array<double, 3> timefactors{
        timefactorGTS, timefactorELTS.result(), timefactorCLTS.result()};
    MPI_Allreduce(
        MPI_IN_PLACE, timefactors.data(), timefactors.size(), MPI_DOUBLE, MPI_SUM, MPI::mpi.comm());

    const auto eltsSpeedup = timefactors[1] / timefactors[0];
    const auto cltsSpeedup = timefactors[2] / timefactors[0];
    logInfo() << "Theoretical speedup to GTS:" << eltsSpeedup << "elementwise LTS;" << cltsSpeedup
              << "clustered LTS (current setup)";

    logInfo() << "Cluster histogram (cell, DR):";
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
