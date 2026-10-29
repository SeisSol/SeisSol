// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_ROCSHMEMKERNELS_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_ROCSHMEMKERNELS_H_

#include <cstddef>
#include <cstdint>
#include <mpi.h>

/**
 * The operations of the SHMEM exchange with ROCSHMEM, as small kernels on a stream that call the
 * device API of ROCSHMEM, as with NVSHMEM (see NvshmemKernels.h): the exchange counts come from
 * device memory when the kernels run, so that a group replayed from a recorded graph counts on and
 * signals and waits for the right values. (The stream API of ROCSHMEM takes these values when an
 * operation gets enqueued, and its puts and waits exist only from version 3.3 on.)
 *
 * The kernels need relocatable device code, linked with the device code of ROCSHMEM in a shared
 * library of their own. All other calls of ROCSHMEM go through it as well: the initialization sets
 * up the state of that device code, and the headers of ROCSHMEM need the HIP compiler.
 */
namespace seissol::solver::rocshmem {

/**
 * Initializes ROCSHMEM on the processes of the communicator; collective. Returns 0 on success.
 */
int initialize(MPI_Comm comm);

void finalize();

/**
 * The processing element of the calling process.
 */
int pe();

/**
 * Allocates symmetric memory; collective.
 */
void* allocate(std::size_t bytes);

/**
 * Frees symmetric memory; collective.
 */
void deallocate(void* pointer);

/**
 * Waits for all processing elements, and completes the operations so far; collective.
 */
void barrierAll();

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
 * The error of the latest kernel launch on the calling thread, and resets it; null if none.
 */
const char* lastError();

} // namespace seissol::solver::rocshmem

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_STREAM_ROCSHMEMKERNELS_H_
