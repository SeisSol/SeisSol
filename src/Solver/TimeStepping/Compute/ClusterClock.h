// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_CLUSTERCLOCK_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_CLUSTERCLOCK_H_

#include "Parallel/Runtime/Stream.h"

namespace seissol::solver {

/**
 * The time at the start of the current step of a cluster, kept where the work on its stream reads
 * it when it runs: kernels read the device copy, host functions the host copy. The clock is set and
 * advanced in the order of the stream, by the same additions that advance the time of the cluster
 * on the host; both thus agree bitwise. Work that is replayed later (e.g. from a graph) reads the
 * time of the step it runs in, not the time it was enqueued with.
 *
 * On a device, one kernel sets or advances both copies; it writes the host copy itself, so that the
 * clock does not wait for the copy engines (e.g. behind a large copy of another stream).
 *
 * Without a device, the clock is a plain value, set and advanced right away.
 */
class ClusterClock {
  public:
  ClusterClock();
  ~ClusterClock();

  ClusterClock(const ClusterClock&) = delete;
  ClusterClock(ClusterClock&&) = delete;
  ClusterClock& operator=(const ClusterClock&) = delete;
  ClusterClock& operator=(ClusterClock&&) = delete;

  /**
   * Sets the time, in the order of the stream.
   */
  void set(double time, parallel::runtime::StreamRuntime& runtime);

  /**
   * Advances the time by the given time step size, in the order of the stream.
   */
  void advance(double timeStepSize, parallel::runtime::StreamRuntime& runtime);

  /**
   * Releases the memory of the clock; needs to happen before the device is finalized.
   */
  void dispose();

  /**
   * The time for host functions on the stream; without a device, the current time.
   */
  [[nodiscard]] const double* host() const { return host_; }

  /**
   * The time for kernels on the stream; null without a device.
   */
  [[nodiscard]] const double* device() const { return device_; }

  private:
  double* host_{nullptr};
  double* device_{nullptr};

#ifdef ACL_DEVICE
  // the address of the host copy on the device
  double* hostOnDevice_{nullptr};
#endif
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_CLUSTERCLOCK_H_
