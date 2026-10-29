// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_EXCHANGEENGINEAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_EXCHANGEENGINEAUX_H_

#include <cstddef>
#include <cstdint>

namespace seissol::kernels::device::aux::exchangeengine {

/// data to copy, in device memory or in host memory that the device maps (by its device address)
struct Region {
  void* target{nullptr};
  const void* source{nullptr};
  std::size_t bytes{0};
};

/**
 * The exchanges of one transport, for the engine (see run()); in device memory.
 *
 * A send gets requested once the copy regions are written; the engine copies them into the host
 * buffer once MPI is done with the send before, and counts that. A receive gets requested once the
 * ghost regions may change; the engine copies the host buffer into them once MPI has received, and
 * counts that.
 */
struct Transport {
  // from the copy regions into the host buffer, and from the host buffer into the ghost regions
  const Region* out{nullptr};
  std::uint32_t outCount{0};
  const Region* in{nullptr};
  std::uint32_t inCount{0};
  // in host memory (by their device addresses): what the engine has copied out and in, and what
  // MPI has sent and received
  std::uint32_t* copiedOut{nullptr};
  const std::uint32_t* sent{nullptr};
  const std::uint32_t* received{nullptr};
  std::uint32_t* copiedIn{nullptr};
  // the sends and receives requested so far (see request()), and those the engine has taken up
  std::uint32_t sendRequested{0};
  std::uint32_t receiveRequested{0};
  std::uint32_t sendClaimed{0};
  std::uint32_t receiveClaimed{0};
};

/**
 * Requests the sends up to `send` of `sender`, and the receives up to `receive` of `receiver`, in
 * the order of the stream; null transports request nothing.
 */
void request(Transport* sender,
             std::uint32_t send,
             Transport* receiver,
             std::uint32_t receive,
             void* stream);

/**
 * Runs the engine for the `count` transports with `blocks` blocks: each one takes up whatever send
 * or receive is ready, in any order, and pauses longer and longer while there is none. It stops
 * once `control[0]` (in host memory, by its device address) is set, and all sends and receives up
 * to `control[1 + t]` and `control[1 + count + t]` of each transport `t` are done. With CUDA, on
 * the stream; with HIP, on an HSA queue of its own (a kernel that keeps running holds up the
 * streams that share its queue, and CLR puts all streams on a few queues), which the
 * synchronization of the device does not wait for. The engine before has to have stopped (see
 * waitForStop()).
 */
void run(Transport* transports,
         std::uint32_t count,
         const std::uint32_t* control,
         unsigned blocks,
         void* stream);

/// waits until the engine that run() has started (with the same stream) has stopped
void waitForStop(void* stream);

} // namespace seissol::kernels::device::aux::exchangeengine

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_EXCHANGEENGINEAUX_H_
