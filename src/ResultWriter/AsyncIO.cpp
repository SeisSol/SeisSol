// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AsyncIO.h"

#include "Parallel/MPI.h"
#include "async/Dispatcher.h"

#include <mpi.h>

namespace seissol::io {

bool AsyncIO::init() {
  async::Dispatcher::init();

  seissol::MPI::mpi.setComm(commWorld());

  return dispatch();
}

void AsyncIO::finalize() {
  // Call parent class
  async::Dispatcher::finalize();

  // Reset the MPI communicator
  seissol::MPI::mpi.setComm(MPI_COMM_WORLD);
}

} // namespace seissol::io
