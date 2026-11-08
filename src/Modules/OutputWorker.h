// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODULES_OUTPUTWORKER_H_
#define SEISSOL_SRC_MODULES_OUTPUTWORKER_H_

#include <async/as/Pin.h>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>

namespace seissol {

/**
 * A thread that processes the synchronization points of the outputs (see
 * Module::processSyncPoint()) while the main thread goes on with the time stepping. Its tasks run
 * one after the other, in the order they were enqueued in. Their MPI calls go over the
 * communicator of the outputs (see Mpi::outputComm()), and their OpenMP regions get teams of their
 * own.
 */
class OutputWorker {
  public:
  /// The thread runs on `cpus`, with up to `ompThreads` threads in its OpenMP regions.
  OutputWorker(const async::as::CpuMask& cpus, int ompThreads);
  ~OutputWorker();

  OutputWorker(const OutputWorker&) = delete;
  OutputWorker& operator=(const OutputWorker&) = delete;
  OutputWorker(OutputWorker&&) = delete;
  OutputWorker& operator=(OutputWorker&&) = delete;

  void enqueue(std::function<void()> task);

  /// The number of CPUs in `cpus`.
  static int cpuCount(const async::as::CpuMask& cpus);

  /// Waits until all tasks enqueued so far have run.
  void wait();

  private:
  void run(const async::as::CpuMask& cpus, int ompThreads);

  std::mutex mutex_;
  std::condition_variable wakeup_;
  std::condition_variable idle_;
  std::deque<std::function<void()>> tasks_;
  bool busy_{false};
  bool stop_{false};
  // the last member, so that the thread starts once the others are there
  std::thread thread_;
};

} // namespace seissol

#endif // SEISSOL_SRC_MODULES_OUTPUTWORKER_H_
