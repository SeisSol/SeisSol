// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "NeighborLinks.h"

#include "Initializer/BasicTypedefs.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Tree/Colormap.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <mpi.h>
#include <set>
#include <vector>

namespace seissol::initializer::internal {

std::map<int, NeighborLink> measureNeighborLinks(const MeshLayout& layout) {
  std::set<int> neighborSet;
  for (const auto& cluster : layout) {
    for (const auto& region : cluster.regions) {
      neighborSet.insert(region.rank);
    }
  }
  const std::vector<int> neighbors(neighborSet.begin(), neighborSet.end());
  const auto count = static_cast<int>(neighbors.size());

  // a communicator of its own, so that no other message can match these
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Comm_dup(Mpi::mpi.comm(), &comm);

  constexpr int SmallBytes = 8;
  constexpr int LargeBytes = 1 << 20;
  constexpr int Repetitions = 3;
  constexpr int Tag = 0;
  const std::vector<char> send(LargeBytes);
  std::vector<std::vector<char>> receive(neighbors.size(), std::vector<char>(LargeBytes));

  // per neighbor, the time until its message arrived, with all of them sent at once after all
  // ranks have started together; the best of a few runs
  const auto round = [&](int bytes) {
    std::vector<double> best(neighbors.size(), std::numeric_limits<double>::infinity());
    std::vector<MPI_Request> requests(2 * neighbors.size(), MPI_REQUEST_NULL);
    std::vector<int> completed(neighbors.size());
    for (int repetition = 0; repetition < Repetitions; ++repetition) {
      Mpi::mpi.barrier(comm);
      const auto start = MPI_Wtime();
      for (int i = 0; i < count; ++i) {
        MPI_Irecv(receive[i].data(), bytes, MPI_BYTE, neighbors[i], Tag, comm, &requests[i]);
      }
      for (int i = 0; i < count; ++i) {
        MPI_Isend(send.data(), bytes, MPI_BYTE, neighbors[i], Tag, comm, &requests[count + i]);
      }
      for (int pending = count; pending > 0;) {
        int done = 0;
        MPI_Waitsome(count, requests.data(), &done, completed.data(), MPI_STATUSES_IGNORE);
        const auto now = MPI_Wtime() - start;
        for (int k = 0; k < done; ++k) {
          best[completed[k]] = std::min(best[completed[k]], now);
        }
        pending -= done;
      }
      MPI_Waitall(count, requests.data() + count, MPI_STATUSES_IGNORE);
    }
    return best;
  };
  const auto small = round(SmallBytes);
  const auto large = round(LargeBytes);
  MPI_Comm_free(&comm);

  std::map<int, NeighborLink> links;
  for (std::size_t i = 0; i < neighbors.size(); ++i) {
    NeighborLink link;
    link.latency = small[i];
    const auto transfer = large[i] - small[i];
    if (transfer > 0) {
      link.bytesPerSecond = static_cast<double>(LargeBytes - SmallBytes) / transfer;
    }
    links[neighbors[i]] = link;
  }
  return links;
}

std::vector<double> haloTimes(const MeshLayout& layout,
                              const LTSColorMap& colormap,
                              const std::map<int, NeighborLink>& links,
                              const std::vector<std::size_t>& bytesPerCell) {
  // what each copy layer waits for, by the rank it comes from: the ghost regions name the copy
  // layer whose cells need them
  std::vector<std::map<int, std::size_t>> bytes(layout.size());
  for (std::size_t color = 0; color < layout.size(); ++color) {
    if (colormap.argument(color).halo != HaloType::Ghost) {
      continue;
    }
    for (const auto& region : layout[color].regions) {
      if (region.remoteId < bytes.size()) {
        bytes[region.remoteId][region.rank] += region.count * bytesPerCell[color];
      }
    }
  }

  std::vector<double> times(layout.size(), 0.0);
  for (std::size_t color = 0; color < layout.size(); ++color) {
    for (const auto& [rank, size] : bytes[color]) {
      const auto link = links.find(rank);
      if (link != links.end()) {
        times[color] = std::max(times[color], link->second.time(size));
      }
    }
  }
  return times;
}

} // namespace seissol::initializer::internal
