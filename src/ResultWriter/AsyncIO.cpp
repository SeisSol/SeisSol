// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AsyncIO.h"

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <async/Config.h>
#include <async/Dispatcher.h>
#include <mpi.h>
#include <string>
#include <utils/env.h>
#include <utils/logger.h>

namespace seissol::io {

AsyncIO::AsyncIO(utils::Env& env)
    : async::Dispatcher(async::Config(utils::Env{env.prefix() + "ASYNC_", "ASYNC_"})),
      prefix_(env.prefix() + "ASYNC_") {}

bool AsyncIO::initDispatcher(const parallel::Pinning& pinning) {
  // where an executor computes as well, the thread that coordinates the output of its node goes
  // where the computation leaves room
  const auto freeCpus = pinning.getFreeCPUsMask();
  const bool roomLeft = !parallel::Pinning::freeCPUsMaskEmpty(freeCpus);
  if (roomLeft) {
    setAffinity(freeCpus);
  }

  // AUTO writes asynchronously where every rank leaves CPUs free for the threads of the output, as
  // it does for a communication thread
  automatic_ = config().mode() == async::Mode::Auto;
  if (automatic_) {
    int room = roomLeft ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &room, 1, MPI_INT, MPI_MIN, seissol::Mpi::mpi.comm());
    if (room == 0) {
      setMode(async::Mode::Sync);
    }
  }

  async::Dispatcher::init();

  const auto layout = async::Dispatcher::layout();
  // the mode as the log names it, with AUTO where that chose it
  const std::string named = automatic_ ? "ASYNC mode AUTO: " : "ASYNC mode ";
  switch (layout.mode) {
  case async::Mode::Sync:
    logInfo() << ("Asynchronous output: off (" + named + "SYNC).").c_str();
    if (automatic_) {
      logInfo() << "Not every rank leaves CPUs free for the threads of the asynchronous output.";
    }
    break;
  case async::Mode::Thread:
    logInfo() << ("Asynchronous output (" + named + "THREAD): a thread on every rank.").c_str();
    break;
  case async::Mode::Auto:
    // init() has settled it
    break;
  case async::Mode::MPI: {
    logInfo() << ("Asynchronous output (" + named + "NODE):").c_str() << layout.executors
              << "executor rank(s) on" << layout.nodes << "node(s)," << layout.executorsPerNode
              << "per group of" << layout.groupNodes << "node(s), each collecting the output of"
              << layout.minGroupRanks << "to" << layout.maxGroupRanks << "ranks;"
              << (layout.exclusive ? "they do not compute." : "they compute as well.");
    if (!layout.rankOrdered) {
      logError() << "The ranks of a node are not numbered one after another (as with a round-robin"
                 << "placement), so the executors of ASYNC mode NODE would not hold the output in"
                 << "the order of the ranks. Place the ranks block-wise, or use ASYNC mode THREAD.";
    }
    if (layout.exclusive) {
      logError() << "SeisSol sets up its output once the simulation runs, which only executors"
                 << "that compute as well can take on. Unset"
                 << (prefix_ + "EXECUTOR_EXCLUSIVE").c_str() << "or use ASYNC mode THREAD.";
    }
    int noRoom = (isExecutor() && !roomLeft) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &noRoom, 1, MPI_INT, MPI_MAX, seissol::Mpi::mpi.comm());
    if (noRoom != 0) {
      logWarning() << "No CPUs are left free on (some) executor ranks: the threads that coordinate"
                   << "and write the output share the cores of the computation there.";
    }
    break;
  }
  }

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
