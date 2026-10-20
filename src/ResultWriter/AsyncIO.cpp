// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AsyncIO.h"

#include "Parallel/MPI.h"

#include <async/Config.h>
#include <async/Dispatcher.h>
#include <mpi.h>
#include <utils/env.h>

namespace seissol::io {

AsyncIO::AsyncIO(utils::Env& env)
    : async::Dispatcher(async::Config(utils::Env{env.prefix() + "ASYNC_", "ASYNC_"})) {}

bool AsyncIO::initDispatcher() {
  async::Dispatcher::init();

  seissol::Mpi::mpi.setComm(commWorld());

  return dispatch();
}

void AsyncIO::finalizeDispatcher() {
  // Call parent class
  async::Dispatcher::finalize();

  // Reset the MPI communicator
  seissol::Mpi::mpi.setComm(MPI_COMM_WORLD);
}

} // namespace seissol::io
