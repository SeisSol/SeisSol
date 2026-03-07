// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ParallelFS.h"

#include "Common/Filesystem.h"
#include "Parallel/MPI.h"

#include <thread>
#include <utils/logger.h>

namespace seissol {

bool entryExistsGlobally(const filesystem::directory_entry& entry) {
  const auto local = directoryExists(entry) ? 1 : 0;
  const auto total = Mpi::mpi.allreduce(local, MPI_MIN);
  return total == 1;
}

void waitExistence(const filesystem::path& path) {
  std::size_t counter = 0;
  auto entry = filesystem::directory_entry(path);
  while (!entryExistsGlobally(entry)) {
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(10ms);
    entry.refresh();
    ++counter;
    if (counter % 10 == 0) {
      logInfo() << "Still waiting for" << path << "to become globally visible.";
    }
  }
}

void moveEntryGlobally(const filesystem::path& entry, const filesystem::path& newentry) {
  logInfo() << "Rename" << entry << "to" << newentry;
  if (Mpi::mpi.rank() == 0) {
    filesystem::rename(entry, newentry);
  }

  waitExistence(newentry);
}

void createDirectoryGlobally(const filesystem::path& entry) {
  logInfo() << "Create directory" << entry;
  if (Mpi::mpi.rank() == 0) {
    filesystem::create_directory(entry);
  }

  // wait for existence propagation to all
  waitExistence(entry);
}

void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp) {
  std::stringstream fullName;
  fullName << fileName << '.' << fileExtension;
  const auto path = seissol::filesystem::path(fullName.str());
  const auto entry = seissol::filesystem::directory_entry(path);

  // file should already exist for a much longer time
  if (entryExistsGlobally(entry)) {
    const auto actualTimeStamp = timeStamp.value_or(
        []() { return utils::TimeUtils::timeAsString("%Y-%m-%dT%H-%M-%S", time(nullptr)); }());
    std::stringstream backupFileName;
    backupFileName << fileName << ".bak_" << actualTimeStamp << '.' << fileExtension;
    const auto copyPath = seissol::filesystem::path(backupFileName.str());
    moveEntryGlobally(seissol::filesystem::directory_entry(path),
                      seissol::filesystem::directory_entry(copyPath));
  }
}

} // namespace seissol
