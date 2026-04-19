// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Datatype.h"
#include "IO/Instance/Point/Csv.h"
#include "IO/Instance/Point/TableWriter.h"

#include <memory>
#include <string>

namespace seissol::unit_test {
using namespace seissol::io::instance::point;
using namespace seissol::io::datatype;

TEST_CASE("Csv quote simple string" * doctest::test_suite("io")) {
  const Csv csv;
  CHECK(csv.quote("hello") == "\"hello\"");
  CHECK(csv.quote("") == "\"\"");
  CHECK(csv.quote("abc") == "\"abc\"");
}

TEST_CASE("Csv quote escapes internal quotes" * doctest::test_suite("io")) {
  const Csv csv;
  std::ostringstream stream;
  csv.quote(stream, "say \"hi\"");
  // Internal quotes should be doubled: say ""hi""
  CHECK(stream.str() == "\"say \"\"hi\"\"\"");
}

TEST_CASE("Csv quote no special chars" * doctest::test_suite("io")) {
  const Csv csv;
  std::ostringstream stream;
  csv.quote(stream, "plain");
  CHECK(stream.str() == "\"plain\"");
}

TEST_CASE("Csv header" * doctest::test_suite("io")) {
  Csv csv;
  auto f64 = std::make_shared<F64Datatype>();
  csv.addQuantity(TableQuantity{"time", f64});
  csv.addQuantity(TableQuantity{"value", f64});

  std::string hdr = csv.header();
  // Expected: "time";"value"\n
  CHECK(hdr == "\"time\";\"value\"\n");
}

TEST_CASE("Csv header single column" * doctest::test_suite("io")) {
  Csv csv;
  auto f32 = std::make_shared<F32Datatype>();
  csv.addQuantity(TableQuantity{"x", f32});
  CHECK(csv.header() == "\"x\"\n");
}

TEST_CASE("Csv rows with data" * doctest::test_suite("io")) {
  Csv csv;
  auto f64 = std::make_shared<F64Datatype>();
  auto i32 = std::make_shared<IntegerDatatype>(sizeof(long long), true);

  csv.addQuantity(TableQuantity{"time", f64});
  csv.addQuantity(TableQuantity{"step", i32});

  // Add a row: time=1.5, step=10
  const double t1 = 1.5;
  const long long s1 = 10;
  csv.addCell(t1);
  csv.addCell(s1);

  // Add another row: time=2.5, step=20
  const double t2 = 2.5;
  const long long s2 = 20;
  csv.addCell(t2);
  csv.addCell(s2);

  const std::string result = csv.rows();

  // Each row should have two quoted values separated by ;
  // Values are produced by toStringRaw, so for double it's %.16g style
  CHECK(result.find("1.5") != std::string::npos);
  CHECK(result.find("10") != std::string::npos);
  CHECK(result.find("2.5") != std::string::npos);
  CHECK(result.find("20") != std::string::npos);

  // Two rows = two newlines
  std::size_t newlines = 0;
  for (const char c : result) {
    if (c == '\n') {
      ++newlines;
    }
  }
  CHECK(newlines == 2);
}

TEST_CASE("Csv rows empty" * doctest::test_suite("io")) {
  Csv csv;
  auto f64 = std::make_shared<F64Datatype>();
  csv.addQuantity(TableQuantity{"x", f64});

  // No cells added → empty rows
  CHECK(csv.rows().empty());
}

TEST_CASE("Csv resetStorage clears data" * doctest::test_suite("io")) {
  Csv csv;
  auto f64 = std::make_shared<F64Datatype>();
  csv.addQuantity(TableQuantity{"x", f64});

  const double val = 42.0;
  csv.addCell(val);
  CHECK_FALSE(csv.rows().empty());

  csv.resetStorage();
  CHECK(csv.rows().empty());
}

TEST_CASE("TableWriter addQuantity and getRowDatatype" * doctest::test_suite("io")) {
  // Use Csv as a concrete TableWriter subclass
  Csv csv;
  auto f64 = std::make_shared<F64Datatype>();
  auto f32 = std::make_shared<F32Datatype>();

  csv.addQuantity(TableQuantity{"alpha", f64});
  csv.addQuantity(TableQuantity{"beta", f32});

  auto rowType = csv.getRowDatatype();
  // Should be a struct with size = 8 + 4 = 12
  CHECK(rowType->size() == 12);
}

} // namespace seissol::unit_test
