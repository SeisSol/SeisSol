// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ParallelFS.h"

#include "Core/Parallel/MPI.h"

#include <filesystem>
#include <mpi.h>
#include <string>
#include <system_error>
#include <utils/logger.h>

namespace seissol {

namespace {
std::string hostName() {
  // (gethostname is POSIX; on Windows, it is part of Winsock and needs extra initialization)
  std::string name(MPI_MAX_PROCESSOR_NAME, '\0');
  int length = 0;
  if (MPI_Get_processor_name(name.data(), &length) != MPI_SUCCESS) {
    return "unknown-host";
  }
  name.resize(length);
  return name;
}
} // namespace

void createDirectoryGlobally(const std::filesystem::path& path) {
  // one rank creates it, so that the metadata load does not grow with the rank count
  if (Mpi::mpi.rank() == 0) {
    std::error_code error;
    const bool created = std::filesystem::create_directories(path, error);
    if (error) {
      logError() << "Could not create the directory" << path.string() << ":" << error.message();
    }
    if (created) {
      logInfo() << "Created the directory" << path.string();
    }
  }

  // no rank looks for the directory before it exists
  Mpi::mpi.barrier();

  std::error_code error;
  const bool visible = std::filesystem::is_directory(path, error);
  if (!visible) {
    logWarning(true) << "The directory" << path.string() << "is not visible on host" << hostName()
                     << ":" << (error ? error.message() : std::string("not a directory"));
  }

  // completes only once every rank has contributed, i.e. after all warnings above are written
  const auto invisibleRanks = Mpi::mpi.allreduce(visible ? 0 : 1, MPI_SUM);
  if (invisibleRanks > 0) {
    logError() << "The directory" << path.string() << "was created, but" << invisibleRanks << "of"
               << Mpi::mpi.size()
               << "ranks cannot see it. The file system did not make it visible on all nodes "
                  "right away; creating it before the run starts usually avoids this.";
  }
}

} // namespace seissol
