// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Instance/Point/Hdf5Table.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "WriterHarness.t.h"

#include <cstddef>
#include <cstring>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io;
using namespace seissol::io::instance::point;

TableQuantity doubleColumn(const std::string& name) {
  return {name, seissol::io::datatype::inferDatatype<double>()};
}

//! What a point of the first group records, and of the second: a quantity more.
std::vector<TableQuantity> small() { return {doubleColumn("v1"), doubleColumn("v2")}; }
std::vector<TableQuantity> large() {
  return {doubleColumn("v1"), doubleColumn("v2"), doubleColumn("p")};
}

//! The value point @p point records for quantity @p quantity at sample @p sample.
double sampleValue(std::size_t point, std::size_t sample, std::size_t quantity) {
  return 100.0 * static_cast<double>(point) + 10.0 * static_cast<double>(sample) +
         static_cast<double>(quantity);
}

} // namespace

TEST_CASE("IO/Hdf5Table: points are split by what they record" * doctest::test_suite("io")) {
  // three points, two of which record the same quantities
  const std::vector<std::vector<TableQuantity>> pointQuantities{small(), large(), small()};
  Hdf5Table table("receivers", pointQuantities, MPI_COMM_SELF, 4);

  const auto& grouping = table.grouping();
  REQUIRE(grouping.groupCount() == 2);
  CHECK(grouping.group[0] == grouping.group[2]);
  CHECK(grouping.group[0] != grouping.group[1]);

  const auto smallGroup = grouping.group[0];
  const auto largeGroup = grouping.group[1];
  CHECK(table.localPointCount(smallGroup) == 2);
  CHECK(table.localPointCount(largeGroup) == 1);
  CHECK(table.sampleSize(smallGroup) == 2 * sizeof(double));
  CHECK(table.sampleSize(largeGroup) == 3 * sizeof(double));
}

TEST_CASE("IO/Hdf5Table: the samples of a point lie together" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;

  const std::vector<std::vector<TableQuantity>> pointQuantities{small(), small()};
  Hdf5Table table("receivers", pointQuantities, MPI_COMM_SELF, 4);
  const auto group = table.grouping().group[0];
  const auto points = table.localPointCount(group);
  REQUIRE(points == 2);

  // two writes carrying a different number of samples each, which is what a sampling interval
  // that does not divide the output interval produces
  const std::vector<std::size_t> samplesPerWrite{3, 2};
  auto plan = table.makeWriter();

  std::size_t written = 0;
  for (std::size_t step = 0; step < samplesPerWrite.size(); ++step) {
    const auto samples = samplesPerWrite[step];
    auto* storage = table.prepare(group, samples);
    for (std::size_t sample = 0; sample < samples; ++sample) {
      for (std::size_t point = 0; point < points; ++point) {
        const std::vector<double> values{sampleValue(point, written + sample, 0),
                                         sampleValue(point, written + sample, 1)};
        std::memcpy(storage + (sample * points + point) * table.sampleSize(group),
                    values.data(),
                    values.size() * sizeof(double));
      }
    }
    // the other group has no samples this round, but still takes part in the write
    for (std::size_t other = 0; other < table.grouping().groupCount(); ++other) {
      if (other != group) {
        static_cast<void>(table.prepare(other, 0));
      }
    }

    auto write = plan(dir.prefix(), step, static_cast<double>(step));
    unit_test::io::runPlan(write, MPI_COMM_SELF);
    written += samples;
  }

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(dir.prefix() + "-receivers.h5");
  hdf5.openGroup("receivers");

  // the dataset grew by the samples of each write, not by one row per write
  const auto values = hdf5.readData<double>("group" + std::to_string(group));
  REQUIRE(values.size() == written * points * 2);

  for (std::size_t sample = 0; sample < written; ++sample) {
    for (std::size_t point = 0; point < points; ++point) {
      const auto base = (sample * points + point) * 2;
      CHECK(values[base] == doctest::Approx(sampleValue(point, sample, 0)));
      CHECK(values[base + 1] == doctest::Approx(sampleValue(point, sample, 1)));
    }
  }

  // and a row of a dataset can be traced back to the point it belongs to
  const auto index = hdf5.readData<std::uint64_t>("Index");
  REQUIRE(index.size() == pointQuantities.size() * 2);
  for (std::size_t point = 0; point < pointQuantities.size(); ++point) {
    CHECK(index[point * 2] == table.grouping().group[point]);
    CHECK(index[point * 2 + 1] == table.grouping().index[point]);
  }

  hdf5.closeGroup();
  hdf5.closeFile();
}

} // namespace seissol::unit_test
