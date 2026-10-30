// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_SHMEMGROUP_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_SHMEMGROUP_H_

#include <cstdint>

/**
 * A group of the SHMEM exchange as one kernel (see ShmemExchangeScheduler): each workgroup sends
 * to one peer or receives from one, as the protocol does it, and all of them at the same time. The
 * tasks and their regions lie in device memory, built once for each direction; the kernels that
 * take them are those of the libraries (see RocshmemKernels.h).
 */
namespace seissol::solver::shmem {

/// one workgroup of the kernel of a group
struct Task {
  /// 1: receives from `peer`, 0: sends to it
  std::uint32_t receive;
  std::int32_t peer;
  /// the regions of the task: `regionCount` from `firstRegion` on
  std::uint32_t firstRegion;
  std::uint32_t regionCount;
  /// receiving: the clear-to-send signal of the direction, which the task sets at the peer;
  /// sending: the one the peer sets here, which the task waits for
  std::uint64_t* clearToSend;
  /// receiving: the arrival signal here, which the task waits for (the count of the exchange times
  /// the number of regions); sending: the one at the peer, which each put adds one to
  std::uint64_t* arrived;
};

/// receiving: from the staging window into a ghost region; sending: from a copy region into the
/// staging window of the peer (a symmetric address)
struct Region {
  void* destination;
  const void* source;
  std::uint64_t bytes;
};

} // namespace seissol::solver::shmem

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_SHMEMGROUP_H_
