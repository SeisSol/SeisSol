// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Reader/Distribution.h"

#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <numeric>
#include <vector>

namespace seissol::unit_test {

namespace distributortest {
using seissol::io::reader::Distributor;

//! @brief Moves @p source, laid out by @p sourceIds, into the order of @p targetIds.
inline std::vector<std::int64_t> moved(const std::vector<std::size_t>& sourceIds,
                                       const std::vector<std::size_t>& targetIds,
                                       const std::vector<std::int64_t>& source) {
  Distributor distributor(MPI_COMM_SELF);
  distributor.setup(sourceIds, targetIds);
  std::vector<std::int64_t> target(targetIds.size(), 0);
  distributor.distribute(target.data(), source.data()).complete();
  return target;
}
} // namespace distributortest

using namespace distributortest;

TEST_CASE("IO/Distributor: data follows the identifiers it is asked for" *
          doctest::test_suite("io")) {
  // what a reader does with a mesh it did not read in the order it needs it
  const std::vector<std::size_t> sourceIds{10, 11, 12, 13};
  const std::vector<std::int64_t> source{100, 101, 102, 103};

  SUBCASE("The same order returns the same data") {
    CHECK(moved(sourceIds, sourceIds, source) == source);
  }

  SUBCASE("A different order rearranges it") {
    const std::vector<std::size_t> targetIds{13, 10, 12, 11};
    CHECK(moved(sourceIds, targetIds, source) == std::vector<std::int64_t>{103, 100, 102, 101});
  }

  SUBCASE("Asking for fewer leaves the rest behind") {
    const std::vector<std::size_t> targetIds{12, 10};
    CHECK(moved(sourceIds, targetIds, source) == std::vector<std::int64_t>{102, 100});
  }

  SUBCASE("Asking for one twice hands it out twice") {
    const std::vector<std::size_t> targetIds{11, 11, 13};
    CHECK(moved(sourceIds, targetIds, source) == std::vector<std::int64_t>{101, 101, 103});
  }

  SUBCASE("Nothing asked for is nothing returned") { CHECK(moved(sourceIds, {}, source).empty()); }
}

TEST_CASE("IO/Distributor: a transform is applied on the way" * doctest::test_suite("io")) {
  const std::vector<std::size_t> sourceIds{0, 1, 2};
  const std::vector<std::size_t> targetIds{2, 1, 0};
  const std::vector<std::int64_t> source{7, 8, 9};

  Distributor distributor(MPI_COMM_SELF);
  distributor.setup(sourceIds, targetIds);

  std::vector<double> target(targetIds.size(), 0.0);
  distributor
      .distributeTransform<double, std::int64_t>(
          target.data(),
          source.data(),
          [](void* to, const void* from) {
            *static_cast<double*>(to) =
                0.5 * static_cast<double>(*static_cast<const int64_t*>(from));
          })
      .complete();

  CHECK(target == std::vector<double>{4.5, 4.0, 3.5});
}

} // namespace seissol::unit_test
