// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AsyncIO.h"

#include "Parallel/MPI.h"

#include <mpi.h>

#include "async/Dispatcher.h"

namespace seissol::io {

bool AsyncIO::init() {
  async::Dispatcher::init();

#ifdef USE_MPI
  seissol::MPI::mpi.setComm(commWorld());
#endif // USE_MPI

  return dispatch();
}

void AsyncIO::finalize() {
  // Call parent class
  async::Dispatcher::finalize();

#ifdef USE_MPI
  // Reset the MPI communicator
  seissol::MPI::mpi.setComm(MPI_COMM_WORLD);
#endif // USE_MPI
}

} // namespace seissol::io
