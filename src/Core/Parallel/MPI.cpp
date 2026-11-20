// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "MPI.h"

#include <mpi.h>
#include <utils/logger.h>

void seissol::Mpi::init(int& argc, char**& argv) {
  // Note: Strictly speaking, we only require MPI_THREAD_MULTIPLE if using
  // a communication thread and/or async I/O.
  // The safer (and more sane) option is to enable it by default.
  const int required = MPI_THREAD_MULTIPLE;
  int provided = 0;
  MPI_Init_thread(&argc, &argv, required, &provided);

  setComm(MPI_COMM_WORLD);

  // Test this after setComm() to get the correct rank_
  if (provided < required) {
    logError() << utils::nospace << "Provided MPI thread support (" << provided
               << ") is less than required thread support (" << required << ").";
  }
}

// (the class wraps the calls that clang-tidy points to it elsewhere)
// NOLINTBEGIN(bugprone-unsafe-functions)
void seissol::Mpi::setComm(MPI_Comm comm) {
  comm_ = comm;
  if (outputComm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&outputComm_);
  }
  MPI_Comm_dup(comm, &outputComm_);

  MPI_Comm_rank(comm, &rank_);
  MPI_Comm_size(comm, &size_);

  if (sharedMemComm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&sharedMemComm_);
  }
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &sharedMemComm_);
  MPI_Comm_rank(sharedMemComm_, &sharedMemMpiRank_);
  MPI_Comm_size(sharedMemComm_, &sharedMemMpiSize_);
}
// NOLINTEND(bugprone-unsafe-functions)

void seissol::Mpi::finalize() {
  if (sharedMemComm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&sharedMemComm_);
  }
  if (outputComm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&outputComm_);
  }
  MPI_Finalize();
}

seissol::Mpi seissol::Mpi::mpi;
