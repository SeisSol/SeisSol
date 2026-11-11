// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Core.h"
#include "Expr/Interp.h"
#include "Expr/Program.h"
#include "Expr/Record.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {

/// The program of the tests of the records: a select, a group, an output with ties.
const char* const RecordProgram = "#!sderiv 1\nin x, group\nout u = select(lt(x, 0.5), x, 2 * x)\n"
                                  "out g = group * 1\n";

/// The identity of the point `point` of the whole point set.
std::uint64_t identityOf(std::size_t point) { return 1000003ULL * point + 7; }

/// The record of `rank` that evaluates the points `points` of the whole point set, in their order.
expr::Record recordOfRank(const expr::Program& program,
                          const std::vector<std::size_t>& points,
                          int rank,
                          bool identified = true) {
  const std::size_t count = points.size();
  std::vector<double> x(count);
  std::vector<double> group(count);
  std::vector<std::uint64_t> identity(count);
  for (std::size_t i = 0; i < count; ++i) {
    identity[i] = identityOf(points[i]);
    x[i] = std::fmod(0.6180339887 * static_cast<double>(points[i]), 1.0);
    group[i] = static_cast<double>(points[i] % 5);
  }
  std::vector<double> u(count);
  std::vector<double> g(count);
  reader::scripting::DataTable table(count);
  table.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
  table.bindViewConst<double>("group", reader::scripting::Direction::In, group.data());
  table.bindView<double>("u", reader::scripting::Direction::Out, u.data());
  table.bindView<double>("g", reader::scripting::Direction::Out, g.data());
  if (identified) {
    table.setIdentity([&identity](std::size_t p) { return identity[p]; });
  }
  auto binding = expr::Binding::bind(program, table);
  reader::datafield::GridStore store;
  expr::BackendOptions options;
  options.quiet = true;
  const auto kernel = expr::makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);
  const auto comparison = expr::compareWithInterpreter(program, binding, store, &table, nullptr);
  REQUIRE(comparison.has_value());
  return expr::recordOf(program, "interpreter", *comparison, {"build = the test"}, rank);
}

std::set<std::uint64_t> identitiesOf(const expr::Record& record) {
  std::set<std::uint64_t> identities;
  for (const auto& point : record.sample) {
    identities.insert(point.identity);
  }
  return identities;
}

/// `first` to `last` (exclusive), in a fixed order other than theirs.
std::vector<std::size_t> shuffled(std::size_t first, std::size_t last) {
  std::vector<std::size_t> points;
  for (std::size_t point = first; point < last; ++point) {
    points.push_back(point);
  }
  std::sort(points.begin(), points.end(), [](std::size_t a, std::size_t b) {
    return expr::splitMix64(a ^ 0x5555U) < expr::splitMix64(b ^ 0x5555U);
  });
  return points;
}

} // namespace

TEST_SUITE("expr records") {
  TEST_CASE("a record prints, reads back and evaluates again to its bits") {
    const expr::Program program = expr::loadCore(RecordProgram, "");
    const expr::Record record = recordOfRank(program, shuffled(0, 300), 3);
    CHECK(record.rank == 3);
    CHECK(record.scope.callPoints == 300);
    CHECK(record.scope.globalIdentities);
    CHECK(record.scope.closedAll);
    CHECK(record.scope.scanned);
    const std::string text = expr::printRecord(record);
    CAPTURE(text);
    CHECK(text.find("identities = global\npool = all points\nscan = all points\n") !=
          std::string::npos);
    const expr::Record read = expr::parseRecord(text, "r3.txt");
    CHECK(expr::printRecord(read) == text);
    CHECK(read.sample.size() == record.sample.size());
    CHECK(read.inputs == std::vector<std::string>{"x", "group"});
    CHECK(read.outputs == std::vector<std::string>{"u", "g"});

    reader::datafield::GridStore store;
    const expr::Replay replay = expr::replayRecord(read, store);
    CHECK(replay.points == record.sample.size());
    CHECK(replay.differing == 0);

    // a value the interpreter does not give
    expr::Record changed = read;
    changed.sample[0].interpreter[0] = std::nextafter(changed.sample[0].interpreter[0], 3.0);
    const expr::Replay again = expr::replayRecord(changed, store);
    CHECK(again.differing == 1);
    REQUIRE(again.differences.size() == 1);
    CHECK(again.differences[0].find("`u` at the point of identity") != std::string::npos);

    // a record whose program is not the one of its fingerprint
    expr::Record other = read;
    other.fingerprint = "0000000000000000";
    CHECK_THROWS_AS(static_cast<void>(expr::programOf(other)), std::runtime_error);
    CHECK_THROWS_AS(static_cast<void>(expr::parseRecord("backend = x\n", "")), std::runtime_error);
  }

  TEST_CASE(
      "the samples of the ranks merge into the sample of all points, whatever the partition") {
    const expr::Program program = expr::loadCore(RecordProgram, "");
    constexpr std::size_t Points = 10000;
    // one rank
    const expr::Record whole = recordOfRank(program, shuffled(0, Points), 0);
    CHECK_FALSE(whole.scope.closedAll);
    CHECK(whole.scope.closed == 4096);
    // three, of which one has more points than its pool, and two evaluate ghosts of the next
    std::vector<std::size_t> first = shuffled(0, 1500);
    std::vector<std::size_t> second = shuffled(1500, 4000);
    const std::vector<std::size_t> third = shuffled(4000, Points);
    for (std::size_t ghost = 1500; ghost < 1600; ++ghost) {
      first.push_back(ghost);
    }
    for (std::size_t ghost = 4000; ghost < 4100; ++ghost) {
      second.push_back(ghost);
    }
    const std::vector<expr::Record> ranks = {recordOfRank(program, first, 0),
                                             recordOfRank(program, second, 1),
                                             recordOfRank(program, third, 2)};
    CHECK(ranks[0].scope.closedAll);
    CHECK_FALSE(ranks[2].scope.closedAll);

    reader::datafield::GridStore store;
    const expr::Record merged = expr::mergeRecords(ranks, store);
    const expr::Record mergedWhole = expr::mergeRecords({whole}, store);
    CHECK(merged.ranks == std::vector<int>{0, 1, 2});
    CHECK(identitiesOf(merged) == identitiesOf(whole));
    CHECK(identitiesOf(mergedWhole) == identitiesOf(whole));
    const auto head = [](const expr::Record& record, const std::string& key) {
      for (const auto& [name, value] : record.head) {
        if (name == key) {
          return value;
        }
      }
      return std::string();
    };
    CHECK(head(merged, "selects") == "exact");
    CHECK(head(merged, "inconsistent").empty());

    // a merged record prints and reads back as well
    const std::string text = expr::printRecord(merged);
    CAPTURE(text);
    CHECK(text.find("merged = 0 1 2\n") != std::string::npos);
    CHECK(expr::printRecord(expr::parseRecord(text, "")) == text);

    // a point with another value on another rank is found
    std::vector<expr::Record> broken = ranks;
    expr::Record::Point copy = ranks[1].sample.front();
    copy.rank = 0;
    copy.interpreter[1] += 1.0;
    broken[0].sample.push_back(copy);
    CHECK(head(expr::mergeRecords(broken, store), "inconsistent").find("1 points") == 0);

    // the records of a consumer without identities do not merge
    const std::vector<expr::Record> local = {recordOfRank(program, first, 0, false),
                                             recordOfRank(program, second, 1, false)};
    CHECK_FALSE(local[0].scope.globalIdentities);
    CHECK_THROWS_AS(static_cast<void>(expr::mergeRecords(local, store)), std::runtime_error);
  }
}

} // namespace seissol::unit_test
