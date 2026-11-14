// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_HOSTCOPYAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_HOSTCOPYAUX_H_

#include <cstddef>

namespace seissol::kernels::device::aux::hostcopy {

/**
 * Copies `bytes` from device memory to host memory that the device maps (`hostOnDevice` is its
 * address on the device), in the order of the stream, with a kernel: unlike a copy, it does not
 * wait for the copy engines, e.g. behind a large copy of another stream.
 */
void copyToHost(void* hostOnDevice, const void* device, std::size_t bytes, void* stream);

} // namespace seissol::kernels::device::aux::hostcopy

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_HOSTCOPYAUX_H_
