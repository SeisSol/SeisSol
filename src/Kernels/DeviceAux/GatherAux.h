// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_GATHERAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_GATHERAUX_H_

#include "Parallel/BlockSet.h"

#include <cstddef>

namespace seissol::kernels::device::aux::gather {

/**
 * Copies blocks of several sizes between their places in device memory and one packed buffer that
 * the device can access, in the order of the stream, with one kernel (for up to 32 sets): block `i`
 * of a set is at `blocks[set.begin + i]` (a null block gets skipped), and at
 * `packed + set.offset + i * set.size` in the packed buffer. Gathers into the packed buffer, or
 * scatters from it (`scatter`). `blocks` is a table in device memory, `sets` one on the host.
 */
void copyBlocks(const parallel::BlockSet* sets,
                std::size_t setCount,
                void* const* blocks,
                void* packed,
                bool scatter,
                void* stream);

} // namespace seissol::kernels::device::aux::gather

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_GATHERAUX_H_
