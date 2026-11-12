// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/GridFactory.h"
#include "Reader/Datafield/ScecGrid.h"
#include "Reader/Datafield/TableGrid.h"

#include <cstdio>
#include <exception>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace seissol::unit_test {

namespace {

namespace df = reader::datafield;

/// The message of what `body` throws; empty where it throws nothing.
template <typename Body>
std::string messageOf(const Body& body) {
  try {
    body();
  } catch (const std::exception& error) {
    return error.what();
  }
  return {};
}

} // namespace

TEST_CASE("Datafield table: the rows of a table are the nodes of its grid, in any order") {
  const auto table = df::parseTextTable("# rho and vs over x and z\n"
                                        "z, x, rho  vs note\n"
                                        "\n"
                                        "10, 2, 0x1.8p+7, -3, 1   # a hexadecimal number\n"
                                        "0 0 100 -0 2\n"
                                        "10 0 110 -1 3\n"
                                        "0 1 200 -1 4\n"
                                        "0, 2, 300, -2, 5\n"
                                        "10 1 210 -2 6\n",
                                        "t.csv");
  CHECK(table.names == std::vector<std::string>{"z", "x", "rho", "vs", "note"});
  REQUIRE(table.rows() == 6);
  CHECK(table.lines.front() == 4);
  CHECK(table.lines.back() == 9);
  // the components in the order of the declaration, the column `note` not read
  const auto data = df::tableGridData(table, {"x", "z"}, {"vs", "rho"}, "t.csv");
  REQUIRE(data.knots.size() == 2);
  CHECK(data.knots[0] == std::vector<double>{0, 1, 2});
  CHECK(data.knots[1] == std::vector<double>{0, 10});
  // the components of a node together, axis 0 (x) fastest
  CHECK(data.values == std::vector<double>{-0.0, 100, -1, 200, -2, 300, -1, 110, -2, 210, -3, 192});
  // of one axis, a profile
  const auto profile = df::tableGridData(df::parseTextTable("depth rho\n"
                                                            "3000 2900\n"
                                                            "0 2600\n"
                                                            "1000 2700\n",
                                                            "p.csv"),
                                         {"depth"},
                                         {"rho"},
                                         "p.csv");
  CHECK(profile.knots == std::vector<std::vector<double>>{{0, 1000, 3000}});
  CHECK(profile.values == std::vector<double>{2600, 2700, 2900});
}

TEST_CASE("Datafield table: a table that is no grid names where") {
  const auto gridOf = [](const std::string& text,
                         const std::vector<std::string>& coordinates = {"x", "z"}) {
    return messageOf([&] {
      static_cast<void>(
          df::tableGridData(df::parseTextTable(text, "t.csv"), coordinates, {"v"}, "t.csv"));
    });
  };
  const std::string head = "x z v\n";
  CHECK(gridOf(head + "0 0 1\n1 0 2\n0 1 3\n0 0 4\n1 1 5\n") ==
        "t.csv:5: the node (x = 0, z = 0) is given twice, before in line 2");
  CHECK(gridOf(head + "0 0 1\n1 0 2\n0 1 3\n") ==
        "the table t.csv gives no node (x = 1, z = 1): its rows are the nodes of a grid, each "
        "once");
  CHECK(gridOf(head + "0 0 1\n", {"x", "y"}) ==
        "the table t.csv has no column `y`; it has x, z, v");
  CHECK(gridOf(head + "0 nan 1\n") == "t.csv:2: the coordinate `z` is nan, not finite");
  CHECK(gridOf(head + "0 0\n") == "t.csv:2: 2 values for 3 columns");
  CHECK(gridOf(head + "0 0 one\n") == "t.csv:2: `one` is no number");
  CHECK(gridOf("x x v\n0 0 1\n") == "t.csv:1: the column `x` is named twice");
  CHECK(gridOf("# nothing\n") == "t.csv: no line names the columns of the table");
  CHECK(gridOf(head) == "the table t.csv has no rows");
  CHECK(messageOf([] { static_cast<void>(df::readTextTable("no-such-table.csv")); }) ==
        "Could not open no-such-table.csv for reading.");
}

TEST_CASE("Datafield table: a grid of a table names its coordinates and its components") {
  df::GridSpec spec;
  spec.kind = "table";
  spec.path = "profile.csv";
  spec.coordinates = {"depth"};
  spec.components = {"rho", "vs"};
  spec.interpolation = {"linear"};
  spec.boundary = {"clamp"};
  const auto desc = df::gridDescOf(spec);
  CHECK(desc.kind == df::GridKind::Table);
  CHECK(desc.columns == spec.components);
  CHECK(desc.coordinates == spec.coordinates);
  CHECK(desc.axisOf(0) == df::AxisKind::Knots);
  // the columns select what the store reads: another selection is another grid
  auto other = spec;
  other.components = {"vs", "rho"};
  CHECK_FALSE(df::gridDescOf(other) == desc);
  CHECK(df::gridDescOf(other).sourceKey() != desc.sourceKey());
  CHECK(df::gridDescOf(other).canonicalKey() == desc.canonicalKey());

  auto withVariable = spec;
  withVariable.variable = "data";
  CHECK(messageOf([&] { static_cast<void>(df::gridDescOf(withVariable)); }) ==
        "a grid of the kind \"table\" has no variable");
  auto withoutCoordinates = spec;
  withoutCoordinates.coordinates.clear();
  CHECK(messageOf([&] { static_cast<void>(df::gridDescOf(withoutCoordinates)); }) ==
        "a grid of the kind \"table\" gives its coordinates");
}

TEST_CASE("Datafield SCEC: a grid has the components a declaration names, in its order") {
  const std::string path = "seissol-test-scec-" + std::to_string(::getpid()) + ".txt";
  {
    // 1 by 1 cells over 10 by 20: four nodes, value c of node (i, j) 100 c + 10 j + i
    std::ofstream file(path);
    file << "1 2\n1 1 10.0 20.0\n0 0 0.0 0.0 1.0 1.0\n";
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        file << i << " " << j << " 0 0";
        for (int c = 0; c < 10; ++c) {
          file << " " << 100 * c + 10 * j + i;
        }
        file << "\n";
      }
    }
  }
  df::GridSpec spec;
  spec.kind = "scec";
  spec.path = path;
  spec.components = {"d_c", "s_normal"};
  spec.interpolation = {"nearest"};
  spec.boundary = {"clamp"};
  const auto desc = df::gridDescOf(spec);
  CHECK(desc.columns == spec.components);
  CHECK(desc.variable.empty());
  const auto grid = df::makeGrid(desc);
  std::remove(path.c_str());
  REQUIRE(grid != nullptr);
  REQUIRE(grid->components() == 2);
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      const double x[2] = {10.0 * i, 20.0 * j};
      double out[2] = {-1.0, -1.0};
      grid->sample(x, out);
      CHECK(out[0] == 700.0 + 10.0 * j + i);
      CHECK(out[1] == 10.0 * j + i);
    }
  }
  auto unknown = spec;
  unknown.components = {"friction"};
  CHECK(messageOf([&] {
          static_cast<void>(df::gridDescOf(unknown));
        }).find("a SCEC file has no component \"friction\"; it has s_normal, s_strike, s_dip") ==
        0);
  auto withVariable = spec;
  withVariable.variable = "z";
  CHECK(messageOf([&] { static_cast<void>(df::gridDescOf(withVariable)); }) ==
        "a grid of the kind \"scec\" has no variable");
}

TEST_CASE("Datafield grid: the values of an inline grid nest as its axes") {
  // values of two axes (3 and 2 knots) and two components, spelled as a source spells them
  const auto number = [](double value) {
    df::SpelledValue spelled;
    spelled.number = value;
    return spelled;
  };
  const auto tuple = [](std::vector<df::SpelledValue> items, std::size_t mark = 0) {
    df::SpelledValue spelled;
    spelled.tuple = true;
    spelled.items = std::move(items);
    spelled.mark = mark;
    return spelled;
  };
  const auto node = [&](double a, double b) { return tuple({number(a), number(b)}); };
  const auto values = tuple({tuple({node(1, 2), node(3, 4), node(5, 6)}),
                             tuple({node(7, 8), node(9, 10), node(11, 12)}, 7)});
  const std::vector<std::vector<double>> knots = {{0, 1, 2}, {10, 20}};
  std::vector<double> flat;
  CHECK_FALSE(df::flattenGridValues(values, knots, {"a", "b"}, flat).has_value());
  CHECK(flat == std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
  // a row of two nodes, which the axis has three of: where it stands
  auto shortRow = values;
  shortRow.items[1].items.pop_back();
  const auto problem = df::flattenGridValues(shortRow, knots, {"a", "b"}, flat);
  REQUIRE(problem.has_value());
  CHECK(problem->mark == 7);
  CHECK(problem->what == "expected the 3 nodes along axis 0, not 2");
}

} // namespace seissol::unit_test
