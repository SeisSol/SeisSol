// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_HOSTCOPYAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_HOSTCOPYAUX_H_

#include <cstddef>
#include <cstdint>

namespace seissol::kernels::device::aux::hostcopy {

/**
 * Copies `bytes` from device memory to host memory that the device maps (`hostOnDevice` is its
 * address on the device), in the order of the stream, with a kernel: unlike a copy, it does not
 * wait for the copy engines, e.g. behind a large copy of another stream.
 */
void copyToHost(void* hostOnDevice, const void* device, std::size_t bytes, void* stream);

/**
 * Copies `bytes` from host memory that the device maps (`hostOnDevice` is its address on the
 * device) to device memory, in the order of the stream, with a kernel (see copyToHost()).
 */
void copyFromHost(void* device, const void* hostOnDevice, std::size_t bytes, void* stream);

/**
 * Writes `value` to a counter in host memory that the device maps (`counterOnDevice` is its
 * address on the device), in the order of the stream: once the host sees the value, it sees all
 * that the stream has written before as well (e.g. with copyToHost()).
 */
void signalHost(std::uint32_t* counterOnDevice, std::uint32_t value, void* stream);

/**
 * Lets the stream wait until a counter in host memory that the device maps (`counterOnDevice` is
 * its address on the device) has reached at least `value`, with a kernel that polls it, pausing
 * longer and longer in between: unlike the wait of HIP (a kernel that polls after the shortest
 * pause), several of them do not flood the memory of the host with reads.
 */
void waitHost(const std::uint32_t* counterOnDevice, std::uint32_t value, void* stream);

} // namespace seissol::kernels::device::aux::hostcopy

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_HOSTCOPYAUX_H_
