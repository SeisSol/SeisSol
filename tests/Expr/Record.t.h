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
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
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

  TEST_CASE("a program with state is compared and recorded from its states before the call") {
    const expr::Program program = expr::loadCore(
        "#!sderiv 1\nin x\nstate s = 0.5\ndef t0 = s + x\nnext s = t0\ndef t1 = 2 * s\n"
        "out u = t1\n",
        "");
    constexpr std::size_t Points = 50;
    std::vector<double> x(Points);
    std::vector<double> u(Points);
    for (std::size_t i = 0; i < Points; ++i) {
      x[i] = 0.25 * static_cast<double>(i);
    }
    reader::scripting::DataTable table(Points);
    table.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
    table.bindView<double>("u", reader::scripting::Direction::Out, u.data());
    auto binding = expr::Binding::bind(program, table);
    reader::datafield::GridStore store;
    expr::BackendOptions options;
    options.quiet = true;
    const auto kernel = expr::makeKernel(program, binding, store, options);
    kernel->precompute(table);
    // the second call, from the states the first one left
    kernel->run(table);
    std::vector<double> before(Points);
    for (std::size_t i = 0; i < Points; ++i) {
      before[i] = 0.5 + x[i];
    }
    kernel->run(table);
    const auto comparison =
        expr::compareWithInterpreter(program, binding, store, &table, nullptr, false, &before);
    REQUIRE(comparison.has_value());
    CHECK(comparison->deviating == 0);
    // without the states before, nothing
    CHECK_FALSE(expr::compareWithInterpreter(program, binding, store, &table, nullptr).has_value());

    const expr::Record record =
        expr::recordOf(program, "interpreter", *comparison, {"build = the test"}, 0);
    CHECK(record.states == std::vector<std::string>{"s"});
    const std::string text = expr::printRecord(record);
    CAPTURE(text);
    CHECK(text.find("point identity x s:state u:interpreter next(s):interpreter\n") !=
          std::string::npos);
    const expr::Record read = expr::parseRecord(text, "");
    CHECK(expr::printRecord(read) == text);
    // the output from the state before, the next state after it: u = 2 s, next s = s + x
    for (const auto& point : read.sample) {
      const double xValue = point.inputs[0];
      const double state = point.inputs[1];
      CHECK(state == 0.5 + xValue);
      CHECK(point.interpreter[0] == 2.0 * state);
      CHECK(point.interpreter[1] == state + xValue);
    }
    const expr::Replay replay = expr::replayRecord(read, store);
    CHECK(replay.points == read.sample.size());
    CHECK(replay.differing == 0);
    // a next state the interpreter does not give
    expr::Record changed = read;
    changed.sample[0].interpreter[1] += 1.0;
    const expr::Replay again = expr::replayRecord(changed, store);
    REQUIRE(again.differing == 1);
    CHECK(again.differences[0].find("`next(s)`") == 0);
  }

  TEST_CASE("the record of a program called repeatedly ends with its last call") {
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() /
        ("seissol-test-last-call-" + std::to_string(std::random_device{}()));
    expr::setRecordDirectory(directory.string(), {"build = the test"});
    const expr::Program program = expr::loadCore(
        "#!sderiv 1\nin x\nstate s = 0.75\ndef t0 = s + x\nnext s = t0\ndef t1 = 3 * s\n"
        "out u = t1\n",
        "");
    constexpr std::size_t Points = 40;
    std::vector<double> x(Points);
    std::vector<double> u(Points);
    for (std::size_t i = 0; i < Points; ++i) {
      x[i] = 0.5 * static_cast<double>(i);
    }
    reader::scripting::DataTable table(Points);
    table.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
    table.bindView<double>("u", reader::scripting::Direction::Out, u.data());
    auto binding = expr::Binding::bind(program, table);
    reader::datafield::GridStore store;
    expr::BackendOptions options;
    options.quiet = true;
    const auto kernel = expr::makeKernel(program, binding, store, options);
    kernel->precompute(table);
    for (int call = 0; call < 3; ++call) {
      kernel->run(table);
    }
    expr::finishRecords();
    expr::setRecordDirectory("");

    std::vector<expr::Record> records;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
      std::ifstream file(entry.path());
      records.push_back(expr::parseRecord(
          std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()),
          entry.path().string()));
    }
    std::filesystem::remove_all(directory);
    REQUIRE(records.size() == 1);
    const expr::Record& record = records[0];
    REQUIRE(record.lastCall.size() == record.sample.size());
    // the third call: from the state the second left, s = 0.75 + 2 x
    for (const auto& point : record.lastCall) {
      const double xValue = point.inputs[0];
      CHECK(point.inputs[1] == 0.75 + 2.0 * xValue);
      CHECK(point.interpreter[0] == 3.0 * (0.75 + 2.0 * xValue));
      CHECK(point.interpreter[1] == 0.75 + 3.0 * xValue);
    }
    bool told = false;
    for (const auto& [key, value] : record.head) {
      told = told || (key == "last call" && value.find("after 2 calls") != std::string::npos);
    }
    CHECK(told);
    const expr::Replay replay = expr::replayRecord(record, store);
    CHECK(replay.points == 2 * record.sample.size());
    CHECK(replay.differing == 0);
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
