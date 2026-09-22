// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Point/Csv.h"
#include "WriterHarness.t.h"

#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io::instance::point;

//! A table of the shape the metadata writers produce: a name, a rank, and a number.
Csv makeTable(CsvFormat format = {}) {
  Csv csv("threads", format);
  csv.addTextColumn("hostname", 16);
  csv.addColumn<std::int32_t>("rank");
  csv.addColumn<double>("weight");
  return csv;
}
} // namespace

TEST_CASE("IO/Csv: only text is quoted by default" * doctest::test_suite("io")) {
  auto csv = makeTable();
  csv.addText("node01");
  csv.addCell<std::int32_t>(3);
  csv.addCell<double>(0.5);

  CHECK(csv.header() == "\"hostname\",\"rank\",\"weight\"\n");

  const auto table = parseCsv(csv.header() + csv.rows());
  REQUIRE(table.rows.size() == 1);
  CHECK(table.header == std::vector<std::string>{"hostname", "rank", "weight"});
  CHECK(table.rows[0][table.column("hostname")] == "node01");
  CHECK(table.rows[0][table.column("rank")] == "3");
}

TEST_CASE("IO/Csv: the punctuation is the writer's to choose" * doctest::test_suite("io")) {
  CsvFormat format;
  format.delimiter = ';';
  format.quoting = CsvQuoting::All;

  auto csv = makeTable(format);
  csv.addText("node01");
  csv.addCell<std::int32_t>(3);
  csv.addCell<double>(0.5);

  const auto text = csv.header() + csv.rows();
  CHECK(text.find(';') != std::string::npos);
  CHECK(text.find("\"3\"") != std::string::npos);

  // and the same settings read it back
  const auto table = parseCsv(text, format);
  REQUIRE(table.rows.size() == 1);
  CHECK(table.rows[0][table.column("rank")] == "3");
}

TEST_CASE("IO/Csv: a value may hold the punctuation" * doctest::test_suite("io")) {
  Csv csv("quoted");
  csv.addTextColumn("text", 24);
  csv.addColumn<std::int32_t>("n");

  csv.addText("a,b");
  csv.addCell<std::int32_t>(1);
  csv.addText("say \"hi\"");
  csv.addCell<std::int32_t>(2);

  const auto table = parseCsv(csv.header() + csv.rows());
  REQUIRE(table.rows.size() == 2);
  CHECK(table.rows[0][0] == "a,b");
  CHECK(table.rows[1][0] == "say \"hi\"");
  CHECK(table.rows[1][1] == "2");
}

TEST_CASE("IO/Csv: text is cut to the length its column holds" * doctest::test_suite("io")) {
  Csv csv("short");
  csv.addTextColumn("text", 4);

  csv.addText("abcdefgh");

  const auto table = parseCsv(csv.header() + csv.rows());
  REQUIRE(table.rows.size() == 1);
  CHECK(table.rows[0][0] == "abcd");
}

TEST_CASE("IO/Csv: a row without a closing newline still counts" * doctest::test_suite("io")) {
  const auto table = parseCsv("a,b\n1,2\n3,4");
  REQUIRE(table.rows.size() == 2);
  CHECK(table.rows[1][1] == "4");
}

TEST_CASE("IO/Csv: a growing table is appended to its file" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  const auto path = dir.prefix() + "-energy.csv";
  {
    // what an earlier run left
    std::ofstream stale(path);
    stale << "stale\n";
  }

  Csv csv("energy");
  csv.addColumn<double>("time");
  csv.addTextColumn("variable", 16);
  csv.addColumn<double>("measurement");

  csv.addCell<double>(0.0);
  csv.addText("elastic_energy");
  csv.addCell<double>(0.1);
  csv.appendFile(path);

  csv.addCell<double>(0.5);
  csv.addText("elastic_energy");
  csv.addCell<double>(1.0 / 3.0);
  csv.appendFile(path);

  const std::ifstream stream(path);
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  // the first write of the run starts the file, so the header is there once and the stale line
  // is gone
  const auto table = parseCsv(buffer.str());
  CHECK(table.header == std::vector<std::string>{"time", "variable", "measurement"});
  REQUIRE(table.rows.size() == 2);
  CHECK(table.rows[0][2] == "0.1");
  CHECK(table.rows[1][0] == "0.5");
  // and a value reads back as the double it was
  CHECK(std::stod(table.rows[1][2]) == 1.0 / 3.0);
}

} // namespace seissol::unit_test
