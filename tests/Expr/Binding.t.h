// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

using reader::scripting::DataTable;
using reader::scripting::Direction;

// Input order is the frontend's first-use order, not alphabetical, and
// Opcode::LoadInput addresses the gathered tile by exactly that index. Tests
// that want a particular channel therefore look it up rather than assuming a
// position -- assuming one would make these tests pass for the wrong reason if
// the order ever changed.
std::size_t inputIndex(const Program& program, const std::string& name) {
  for (std::size_t i = 0; i < program.inputs().size(); ++i) {
    if (program.inputs()[i].name == name) {
      return i;
    }
  }
  FAIL("the program does not read " << name);
  return 0;
}

} // namespace

TEST_SUITE("ExprBinding") {

  TEST_CASE("gather converts every column type into the compute tile") {
    const Program program = compileSderiv("x * 2.0 + y + group", "out");
    constexpr std::size_t NumPoints = 5;

    std::vector<float> x = {1.0F, 2.0F, 3.0F, 4.0F, 5.0F};
    std::vector<double> y = {10.0, 20.0, 30.0, 40.0, 50.0};
    std::vector<std::int32_t> group = {7, 7, 7, 7, 7};
    std::vector<double> out(NumPoints, -1.0);

    DataTable table(NumPoints);
    table.bindView<float>("x", Direction::In, x.data());
    table.bindView<double>("y", Direction::In, y.data());
    table.bindView<std::int32_t>("group", Direction::In, group.data());
    table.bindView<double>("out", Direction::Out, out.data());

    const Binding binding = Binding::bind(program, table);
    REQUIRE(binding.numPoints() == NumPoints);
    REQUIRE(binding.inputs().size() == program.inputs().size());

    std::vector<double> tile(program.inputs().size() * NumPoints, 0.0);
    binding.gather(table, 0, NumPoints, tile.data());

    const std::size_t xi = inputIndex(program, "x");
    const std::size_t yi = inputIndex(program, "y");
    const std::size_t gi = inputIndex(program, "group");
    for (std::size_t lane = 0; lane < NumPoints; ++lane) {
      CHECK(tile[xi * NumPoints + lane] == doctest::Approx(x[lane]));
      CHECK(tile[yi * NumPoints + lane] == doctest::Approx(y[lane]));
      CHECK(tile[gi * NumPoints + lane] == doctest::Approx(group[lane]));
    }

    std::vector<double> result(NumPoints);
    for (std::size_t lane = 0; lane < NumPoints; ++lane) {
      result[lane] = 100.0 + static_cast<double>(lane);
    }
    binding.scatter(table, 0, NumPoints, result.data());
    CHECK(out[0] == doctest::Approx(100.0));
    CHECK(out[4] == doctest::Approx(104.0));
  }

  TEST_CASE("the f32 compute path gathers, scatters and allocates state") {
    Program program = compileSderiv("x + 1.0", "out");
    program.setComputeType(ComputeType::F32);
    constexpr std::size_t NumPoints = 4;

    std::vector<double> x = {1.5, 2.5, 3.5, 4.5};
    std::vector<float> out(NumPoints, -1.0F);

    DataTable table(NumPoints);
    table.bindView<double>("x", Direction::In, x.data());
    table.bindView<float>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);

    std::vector<float> tile(program.inputs().size() * NumPoints, 0.0F);
    binding.gather(table, 0, NumPoints, tile.data());
    const std::size_t xi = inputIndex(program, "x");
    CHECK(tile[xi * NumPoints + 0] == 1.5F);
    CHECK(tile[xi * NumPoints + 3] == 4.5F);

    const std::vector<float> result = {9.0F, 8.0F, 7.0F, 6.0F};
    binding.scatter(table, 0, NumPoints, result.data());
    CHECK(out[0] == 9.0F);
    CHECK(out[3] == 6.0F);

    binding.allocatePersistent(program, 3);
    CHECK(binding.persistentSlotCount() == 3);
    REQUIRE(binding.persistentF32() != nullptr);
  }

  TEST_CASE("an unsorted group column is partitioned and permuted") {
    const Program program = compileSderiv("x + group", "out");
    constexpr std::size_t NumPoints = 6;

    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<std::int32_t> group = {2, 1, 2, 1, 3, 1};
    std::vector<double> out(NumPoints, -1.0);

    DataTable table(NumPoints);
    table.bindView<double>("x", Direction::In, x.data());
    table.bindView<std::int32_t>("group", Direction::In, group.data());
    table.bindView<double>("out", Direction::Out, out.data());

    const Binding binding = Binding::bind(program, table);
    REQUIRE(binding.groupRanges().size() == 3);
    REQUIRE(binding.permutation().size() == NumPoints);

    for (const auto& range : binding.groupRanges()) {
      for (std::size_t k = range.begin; k < range.end; ++k) {
        CHECK(group[binding.permutation()[k]] == range.group);
      }
    }

    // The round trip is the property that matters: a kernel working in permuted
    // space must land its results back on the points the consumer handed in.
    const std::size_t xi = inputIndex(program, "x");
    std::vector<double> tile(program.inputs().size() * NumPoints, 0.0);
    binding.gather(table, 0, NumPoints, tile.data());

    std::vector<double> result(NumPoints);
    for (std::size_t lane = 0; lane < NumPoints; ++lane) {
      result[lane] = tile[xi * NumPoints + lane] * 10.0;
    }
    binding.scatter(table, 0, NumPoints, result.data());
    for (std::size_t p = 0; p < NumPoints; ++p) {
      CHECK(out[p] == doctest::Approx(x[p] * 10.0));
    }
  }

  TEST_CASE("an already-grouped point set keeps its order") {
    const Program program = compileSderiv("x + group", "out");
    constexpr std::size_t NumPoints = 4;

    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<std::int32_t> group = {1, 1, 5, 5};
    std::vector<double> out(NumPoints, 0.0);

    DataTable table(NumPoints);
    table.bindView<double>("x", Direction::In, x.data());
    table.bindView<std::int32_t>("group", Direction::In, group.data());
    table.bindView<double>("out", Direction::Out, out.data());

    const Binding binding = Binding::bind(program, table);
    CHECK(binding.groupRanges().size() == 2);
    // An identity permutation is dropped rather than stored: keeping it would
    // cost an indirection per lane in every gather for no reordering.
    CHECK(binding.permutation().empty());
  }

  TEST_CASE("a program without a group input is not partitioned") {
    const Program program = compileSderiv("x + 1.0", "out");
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> out(2, 0.0);

    DataTable table(2);
    table.bindView<double>("x", Direction::In, x.data());
    table.bindView<double>("out", Direction::Out, out.data());

    const Binding binding = Binding::bind(program, table);
    CHECK(binding.groupRanges().empty());
    CHECK(binding.permutation().empty());
  }

  TEST_CASE("bind rejects tables it cannot resolve unambiguously") {
    const Program program = compileSderiv("x + y", "out");
    std::vector<double> x = {1.0};
    std::vector<double> y = {2.0};
    std::vector<double> out = {0.0};

    SUBCASE("a required input has no column") {
      DataTable table(1);
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("out", Direction::Out, out.data());
      CHECK_THROWS_AS(Binding::bind(program, table), std::invalid_argument);
    }

    SUBCASE("a column name appears twice") {
      DataTable table(1);
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("y", Direction::In, y.data());
      table.bindView<double>("out", Direction::Out, out.data());
      CHECK_THROWS_AS(Binding::bind(program, table), std::invalid_argument);
    }

    SUBCASE("an output is bound to an In-only column") {
      DataTable table(1);
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("y", Direction::In, y.data());
      table.bindView<double>("out", Direction::In, out.data());
      CHECK_THROWS_AS(Binding::bind(program, table), std::invalid_argument);
    }

    // Direction and writability are separate questions: a const view can carry
    // InOut and still have no setter, which would be a null std::function call
    // at the first tile rather than a diagnostic.
    SUBCASE("an output is bound to a read-only view") {
      DataTable table(1);
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("y", Direction::In, y.data());
      table.bindViewConst<double>("out", Direction::InOut, out.data());
      CHECK_THROWS_AS(Binding::bind(program, table), std::invalid_argument);
    }

    SUBCASE("the point set is empty") {
      DataTable table(0);
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("y", Direction::In, y.data());
      table.bindView<double>("out", Direction::Out, out.data());
      CHECK_THROWS_AS(Binding::bind(program, table), std::invalid_argument);
    }

    // Explicitly NOT an error: the consumer may offer more than the program
    // reads, and Binding.h says so.
    SUBCASE("an extra column is accepted") {
      DataTable table(1);
      table.bindView<double>("x", Direction::In, x.data());
      table.bindView<double>("y", Direction::In, y.data());
      table.bindView<double>("out", Direction::Out, out.data());
      table.bindView<double>("unused", Direction::In, x.data());
      CHECK_NOTHROW(Binding::bind(program, table));
    }
  }

} // TEST_SUITE

} // namespace seissol::expr::test
