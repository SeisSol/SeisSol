// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/BatchRecorders/DataTypes/ScratchRanges.h"

#include <array>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Scratch ranges are consecutive and disjoint" * doctest::test_suite("initializer")) {
  constexpr std::size_t EntrySize = 5;

  // e.g. the number of analytical (cell, face) pairs per local face;
  // LTS::AnalyticScratch holds their sum (deriveRequiredScratchpadMemoryForWp)
  const std::array<std::size_t, 4> counts{3, 0, 1, 2};
  std::vector<double> scratch(6 * EntrySize);

  recording::ScratchRanges<double> ranges(scratch.data(), EntrySize);

  // the entries of a face directly follow those of the previous faces; hence, no two
  // (face, index) pairs share an entry, and all of them fit into the scratch
  std::size_t expected = 0;
  for (std::size_t face = 0; face < counts.size(); ++face) {
    CAPTURE(face);
    const double* range = ranges.take(counts[face]);
    for (std::size_t i = 0; i < counts[face]; ++i) {
      CAPTURE(i);
      CHECK(range + i * EntrySize == scratch.data() + expected * EntrySize);
      ++expected;
    }
  }
  CHECK(expected * EntrySize == scratch.size());
}

TEST_CASE("Empty scratch ranges take no space" * doctest::test_suite("initializer")) {
  constexpr std::size_t EntrySize = 3;
  std::vector<float> scratch(2 * EntrySize);

  recording::ScratchRanges<float> ranges(scratch.data(), EntrySize);

  CHECK(ranges.take(0) == scratch.data());
  CHECK(ranges.take(1) == scratch.data());
  CHECK(ranges.take(0) == scratch.data() + EntrySize);
  CHECK(ranges.take(1) == scratch.data() + EntrySize);
  CHECK(ranges.take(0) == scratch.data() + scratch.size());
}

} // namespace seissol::unit_test
