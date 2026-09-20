// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Instance/Point/Grouping.h"

#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <set>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io::instance::point;
using seissol::io::datatype::inferDatatype;

TableQuantity real(const std::string& name) { return {name, inferDatatype<double>()}; }

//! what an elastic receiver records
std::vector<TableQuantity> elastic() {
  return {real("v1"), real("v2"), real("v3"), real("xx"), real("yy"), real("zz")};
}

//! the same plus the pressure a poroelastic one has
std::vector<TableQuantity> poroelastic() {
  auto quantities = elastic();
  quantities.push_back(real("p"));
  return quantities;
}

} // namespace

TEST_CASE("IO/Grouping: a quantity set survives its key" * doctest::test_suite("io")) {
  const auto quantities = poroelastic();
  const auto restored = quantitySetFromKey(quantitySetKey(quantities));

  REQUIRE(restored.size() == quantities.size());
  for (std::size_t i = 0; i < quantities.size(); ++i) {
    CHECK(restored[i].name == quantities[i].name);
    CHECK(restored[i].datatype->size() == quantities[i].datatype->size());
  }

  // the key has to separate sets that differ, in name and in type
  CHECK(quantitySetKey(elastic()) != quantitySetKey(poroelastic()));
  CHECK(quantitySetKey({real("a")}) != quantitySetKey({real("b")}));
  CHECK(quantitySetKey({real("a")}) != quantitySetKey({{"a", inferDatatype<std::int64_t>()}}));
}

TEST_CASE("IO/Grouping: points are gathered by their quantity set" * doctest::test_suite("io")) {
  const std::vector<std::vector<TableQuantity>> points{
      elastic(), poroelastic(), elastic(), elastic(), poroelastic()};

  const auto grouping = groupPoints(points, MPI_COMM_SELF);

  REQUIRE(grouping.groupCount() == 2);
  // points of the same set land in the same group, points of different sets do not
  CHECK(grouping.group[0] == grouping.group[2]);
  CHECK(grouping.group[0] == grouping.group[3]);
  CHECK(grouping.group[1] == grouping.group[4]);
  CHECK(grouping.group[0] != grouping.group[1]);

  REQUIRE(grouping.globalCount.size() == 2);
  CHECK(grouping.globalCount[grouping.group[0]] == 3);
  CHECK(grouping.globalCount[grouping.group[1]] == 2);

  // every group declares the quantities its points share
  CHECK(grouping.quantities[grouping.group[0]].size() == elastic().size());
  CHECK(grouping.quantities[grouping.group[1]].size() == poroelastic().size());

  // inside a group the indices are a permutation of 0..count-1
  for (std::size_t group = 0; group < grouping.groupCount(); ++group) {
    std::set<std::size_t> seen;
    for (std::size_t point = 0; point < points.size(); ++point) {
      if (grouping.group[point] == group) {
        seen.emplace(grouping.index[point]);
      }
    }
    REQUIRE(seen.size() == grouping.globalCount[group]);
    CHECK(*seen.begin() == 0);
    CHECK(*seen.rbegin() == grouping.globalCount[group] - 1);
  }
}

TEST_CASE("IO/Grouping: one set means one group" * doctest::test_suite("io")) {
  const std::vector<std::vector<TableQuantity>> points{elastic(), elastic(), elastic()};
  const auto grouping = groupPoints(points, MPI_COMM_SELF);

  REQUIRE(grouping.groupCount() == 1);
  CHECK(grouping.index == std::vector<std::size_t>{0, 1, 2});
  CHECK(grouping.globalCount[0] == 3);
}

TEST_CASE("IO/Grouping: no points at all" * doctest::test_suite("io")) {
  const auto grouping = groupPoints({}, MPI_COMM_SELF);
  CHECK(grouping.groupCount() == 0);
  CHECK(grouping.index.empty());
}

} // namespace seissol::unit_test
