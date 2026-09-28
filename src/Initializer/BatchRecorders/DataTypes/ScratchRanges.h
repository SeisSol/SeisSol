// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_SCRATCHRANGES_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_SCRATCHRANGES_H_

#include <cstddef>

namespace seissol::recording {

/**
 * Splits a scratchpad into consecutive, disjoint ranges of entries (of entrySize elements each),
 * in the order in which they are taken.
 *
 * Meant for scratchpads which a host function fills while a kernel may still read another range,
 * e.g. the per-face ranges of LTS::AnalyticScratch. The recorder and the code that fills the
 * scratchpad agree on the ranges by taking them in the same order.
 */
template <typename T>
class ScratchRanges {
  public:
  ScratchRanges(T* scratch, std::size_t entrySize) : next_(scratch), entrySize_(entrySize) {}

  /**
   * Returns the start of the next range of count entries.
   */
  T* take(std::size_t count) {
    T* range = next_;
    next_ += count * entrySize_;
    return range;
  }

  private:
  T* next_;
  std::size_t entrySize_;
};

} // namespace seissol::recording

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_SCRATCHRANGES_H_
