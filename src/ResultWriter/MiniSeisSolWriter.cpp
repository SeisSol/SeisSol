// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MiniSeisSolWriter.h"

#include "Common/Filesystem.h"
#include "IO/Instance/Point/Csv.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

namespace {
//! @brief How wide a column has to be to hold every one of @p values.
std::size_t widest(const std::vector<std::string>& values) {
  std::size_t width = 1;
  for (const auto& value : values) {
    width = std::max(width, value.size());
  }
  return width;
}
} // namespace

void seissol::writer::MiniSeisSolWriter::write(double elapsedTime, double weight) {
  auto elapsedTimeVector = seissol::Mpi::mpi.collect(elapsedTime);
  auto weightVector = seissol::Mpi::mpi.collect(weight);

  auto localRanks = seissol::Mpi::mpi.collect(seissol::Mpi::mpi.sharedMemMpiRank());

  if (seissol::Mpi::mpi.rank() == 0) {
    std::vector<size_t> ranks(seissol::Mpi::mpi.size());
    for (size_t i = 0; i < ranks.size(); ++i) {
      ranks[i] = i;
    }

    // the slowest rank first, which is the one a reader of this file is after
    std::sort(ranks.begin(), ranks.end(), [&elapsedTimeVector](const size_t& i, const size_t& j) {
      return elapsedTimeVector[i] > elapsedTimeVector[j];
    });

    const auto& hostNames = seissol::Mpi::mpi.getHostNames();

    seissol::io::instance::point::Csv table("miniSeissol");
    table.addTextColumn("hostname", widest(hostNames));
    table.addColumn<std::size_t>("rank");
    table.addColumn<int>("localRank");
    table.addColumn<double>("elapsedTime");
    table.addColumn<double>("weight");

    for (auto rank : ranks) {
      table.addText(hostNames[rank]);
      table.addCell<std::size_t>(rank);
      table.addCell<int>(localRanks[rank]);
      table.addCell<double>(elapsedTimeVector[rank]);
      table.addCell<double>(weightVector[rank]);
    }

    seissol::filesystem::path path(outputDirectory_);
    path += seissol::filesystem::path("-miniSeissol.csv");
    table.writeFile(path.string());
  }
}
