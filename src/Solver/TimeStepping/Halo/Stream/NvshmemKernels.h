// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_NVSHMEMKERNELS_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_NVSHMEMKERNELS_H_

#include "Solver/TimeStepping/Halo/Stream/ShmemGroup.h"

#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <vector>

/**
 * The operations of the SHMEM exchange with NVSHMEM, as small kernels on a stream that call the
 * device API of NVSHMEM. The exchange counts come from device memory when the kernels run: each
 * group advances the count of its direction by a kernel as well, so that a group replayed from a
 * recorded graph counts on and signals and waits for the right values. (The host API of NVSHMEM
 * takes these values when an operation gets enqueued; a graph would keep them.)
 *
 * The kernels need relocatable device code and the static device library of NVSHMEM, and so does
 * the initialization of NVSHMEM, which sets up the state of the device library for them.
 */
namespace seissol::solver::nvshmem {

/**
 * Initializes NVSHMEM on the processes of the communicator, with its device library; collective.
 * Returns 0 on success.
 */
int initialize(MPI_Comm comm);

void finalize();

/**
 * Adds one to the count.
 */
void advance(std::uint64_t* count, void* stream);

/**
 * Sets the signal on the processing element to the count.
 */
void signalCount(std::uint64_t* signal, const std::uint64_t* count, int pe, void* stream);

/**
 * Waits until the (local) signal has reached the count times the factor.
 */
void waitCount(std::uint64_t* signal,
               const std::uint64_t* count,
               std::uint64_t factor,
               void* stream);

/**
 * Puts the data into the symmetric memory of the processing element, and then adds one to the
 * signal there; non-blocking.
 */
void putSignal(void* destination,
               const void* source,
               std::size_t bytes,
               std::uint64_t* signal,
               int pe,
               void* stream);

/**
 * Waits until all puts have completed.
 */
void quiet(void* stream);

/**
 * A whole group as one kernel (see ShmemGroup.h): a block for each of the `taskCount` tasks, so
 * that the waits for the peers overlap; the count of the direction in device memory gets advanced
 * by the kernel (by its last block, with the help of the counter `finished`, which it leaves at 0).
 */
void group(const shmem::Task* tasks,
           std::size_t taskCount,
           const shmem::Region* regions,
           std::uint64_t windowStride,
           std::uint64_t* count,
           std::uint32_t* finished,
           void* stream);

/**
 * The error of the latest kernel launch on the calling thread, and resets it; null if none.
 */
const char* lastError();

/**
 * Registers the memory that puts read from, `count` ranges of `bytes` bytes at `sources`, with
 * NVSHMEM: each allocation that holds a part of a range, once, as a whole. Puts to another node
 * read their source through the network card (IBRC, IBGDA), which needs the memory registered
 * unless it is symmetric; within a node, the devices read each other's memory directly.
 *
 * Fills `registered` with what got registered, for unregisterSources(); returns null on success,
 * or else what failed.
 */
const char* registerSources(const void* const* sources,
                            const std::size_t* bytes,
                            std::size_t count,
                            std::vector<void*>& registered);

/**
 * Takes back what registerSources() registered.
 */
void unregisterSources(const std::vector<void*>& registered);

} // namespace seissol::solver::nvshmem

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_NVSHMEMKERNELS_H_
