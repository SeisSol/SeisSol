// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_HALOCHECK_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_HALOCHECK_H_

#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace seissol::solver {

/**
 * A transport of the halo exchange, with the regions it moves, and what it connects (for the
 * messages).
 */
struct CheckedTransport {
  HaloTransport* transport{nullptr};
  const RemoteClusterPair* regions{nullptr};
  std::string description;
};

/**
 * Where the regions of the halo exchange reside: in host memory, or in device memory, which the
 * check reaches by copies.
 */
enum class HaloMemory : std::uint8_t { Host, Device };

/**
 * What the check of the halo exchange found on one process.
 */
struct HaloCheckResult {
  /// the ghost regions of the process, and the bytes they received
  std::size_t regions{0};
  std::size_t bytes{0};
  /// the ghost regions that did not receive the data of their sender
  std::size_t wrongRegions{0};
  /// what was wrong with the first of them
  std::string firstWrong;
};

/**
 * The value that a copy region holds at a position for the check: it tells the sending process,
 * the tag of the region and the position. It is an integer below 2^23, which single precision
 * holds exactly, and never zero, so that a ghost region that has received nothing does not pass.
 */
[[nodiscard]] double haloCheckValue(int sender, std::size_t tag, std::size_t position);

/**
 * Tests the halo exchange through the transports of the run, before the time stepping: the copy
 * regions get the values of haloCheckValue(), one exchange in each direction moves them, and each
 * ghost region then has to hold the values of its sender. Afterwards, the copy and the ghost
 * regions are zero, as after their allocation; the first time step writes them anew before
 * anything reads them.
 *
 * `rank` is the rank of this process. Collective over the processes that exchange with each other.
 * If the exchange has not completed after `timeout` (none if zero), it reports what is still in
 * flight, and aborts with `abortOnTimeout`; or else, reports again after each further timeout.
 */
HaloCheckResult checkHaloExchange(const std::vector<CheckedTransport>& transports,
                                  int rank,
                                  HaloMemory memory,
                                  std::chrono::seconds timeout,
                                  bool abortOnTimeout);

/**
 * Checks the halo exchange of each process with checkHaloExchange(), prints what the check found
 * on all processes, and aborts if any ghost region received wrong data. `mode` is the name of the
 * transfer mode, for the messages. Collective over all processes; nothing to do without ghost
 * regions.
 */
void testHaloExchange(const std::vector<CheckedTransport>& transports,
                      HaloMemory memory,
                      const std::string& mode,
                      std::chrono::seconds timeout,
                      bool abortOnTimeout);

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_HALOCHECK_H_
