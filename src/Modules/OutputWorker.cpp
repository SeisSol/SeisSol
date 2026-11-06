// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "OutputWorker.h"

#include "Parallel/MPI.h"

#include <algorithm>
#include <async/as/Pin.h>
#include <functional>
#include <mpi.h>
#include <mutex>
#include <pthread.h>
#include <sched.h>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seissol {

OutputWorker::OutputWorker(const async::as::CpuMask& cpus, int ompThreads)
    : thread_([this, cpus, ompThreads] { run(cpus, ompThreads); }) {}

OutputWorker::~OutputWorker() {
  {
    const std::scoped_lock lock(mutex_);
    stop_ = true;
  }
  wakeup_.notify_all();
  thread_.join();
}

void OutputWorker::enqueue(std::function<void()> task) {
  {
    const std::scoped_lock lock(mutex_);
    tasks_.push_back(std::move(task));
  }
  wakeup_.notify_one();
}

int OutputWorker::cpuCount(const async::as::CpuMask& cpus) {
#ifndef __APPLE__
  return CPU_COUNT(&cpus.set);
#else
  return 1;
#endif
}

void OutputWorker::wait() {
  std::unique_lock lock(mutex_);
  idle_.wait(lock, [this] { return tasks_.empty() && !busy_; });
}

void OutputWorker::run(const async::as::CpuMask& cpus, int ompThreads) {
  cpus.setaffinity_np(pthread_self());
#ifdef _OPENMP
  omp_set_num_threads(std::max(1, ompThreads));
#endif
  Mpi::setThreadComm(Mpi::mpi.outputComm());

  std::unique_lock lock(mutex_);
  while (true) {
    wakeup_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
    // the tasks still enqueued run before the thread stops
    if (tasks_.empty()) {
      break;
    }
    auto task = std::move(tasks_.front());
    tasks_.pop_front();
    busy_ = true;
    lock.unlock();
    task();
    lock.lock();
    busy_ = false;
    idle_.notify_all();
  }

  Mpi::setThreadComm(MPI_COMM_NULL);
}

} // namespace seissol
