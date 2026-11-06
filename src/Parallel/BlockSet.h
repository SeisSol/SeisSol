// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_BLOCKSET_H_
#define SEISSOL_SRC_PARALLEL_BLOCKSET_H_

#include <cstddef>

namespace seissol::parallel {

/**
 * Blocks of one size in a packed buffer (see DataCollectorUntyped and
 * kernels::device::aux::gather::copyBlocks()).
 */
struct BlockSet {
  // the first block of the set in the table of blocks
  std::size_t begin{0};
  std::size_t count{0};
  // the bytes of each block
  std::size_t size{0};
  // where the first block of the set lies in the packed buffer, in bytes
  std::size_t offset{0};
};

} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_BLOCKSET_H_
