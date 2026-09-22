// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/Writer/File/Hdf5Writer.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Dimension.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "WriterHarness.t.h"

#include <fstream>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::unit_test::runfiles {
using namespace seissol::io;

//! A plan writing @p values into a dataset of fixed size, as a checkpoint or a snapshot does.
inline writer::Writer fixedPlan(const std::string& path, const std::vector<double>& values) {
  writer::Writer plan;
  plan.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
      writer::instructions::Hdf5Location(path, {"data"}),
      "values",
      writer::WriteInline::createArray<double>({values.size()}, values),
      datatype::inferDatatype<double>()));
  return plan;
}

//! A plan appending @p value to a dataset that grows with every write, as a time series does.
inline writer::Writer appendPlan(const std::string& path, double value) {
  writer::Writer plan;
  plan.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
      writer::instructions::Hdf5Location(path, {"data"}),
      "series",
      writer::WriteInline::createShaped<double>({writer::Dimension::appended(1)}, {value}),
      datatype::inferDatatype<double>()));
  return plan;
}

inline std::vector<double> readBack(const std::string& path, const std::string& name) {
  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(path);
  hdf5.openGroup("data");
  auto values = hdf5.readData<double>(name);
  hdf5.closeGroup();
  hdf5.closeFile();
  return values;
}

TEST_CASE("IO/RunFiles: a new run replaces the files an earlier run left" *
          doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  const auto path = dir.prefix() + "-checkpoint.h5";

  writer::file::RunFiles first;
  auto plan = fixedPlan(path, {1.0, 1.0});
  unit_test::io::runPlan(plan, MPI_COMM_SELF, &first);

  // the dataset is there already, but it belongs to the earlier run
  writer::file::RunFiles second;
  auto again = fixedPlan(path, {2.0, 2.0, 2.0});
  unit_test::io::runPlan(again, MPI_COMM_SELF, &second);

  CHECK(readBack(path, "values") == std::vector<double>{2.0, 2.0, 2.0});
}

TEST_CASE("IO/RunFiles: a run appends to what it wrote itself" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  const auto path = dir.prefix() + "-series.h5";

  writer::file::RunFiles run;
  for (const double value : {1.0, 2.0, 3.0}) {
    auto plan = appendPlan(path, value);
    unit_test::io::runPlan(plan, MPI_COMM_SELF, &run);
  }
  CHECK(readBack(path, "series") == std::vector<double>{1.0, 2.0, 3.0});

  // a run resuming from a checkpoint continues the file of the run before it
  writer::file::RunFiles resumed;
  resumed.resumed = true;
  auto next = appendPlan(path, 4.0);
  unit_test::io::runPlan(next, MPI_COMM_SELF, &resumed);
  CHECK(readBack(path, "series") == std::vector<double>{1.0, 2.0, 3.0, 4.0});

  // and a new run starts it over
  writer::file::RunFiles fresh;
  auto first = appendPlan(path, 5.0);
  unit_test::io::runPlan(first, MPI_COMM_SELF, &fresh);
  CHECK(readBack(path, "series") == std::vector<double>{5.0});
}

TEST_CASE("IO/RunFiles: a binary write that does not append starts the file over" *
          doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  const auto path = dir.prefix() + "-payload.bin";

  const auto write = [&](const std::string& text, bool append) {
    writer::Writer plan;
    plan.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
        path, writer::WriteInline::createString(text), 0, append));
    unit_test::io::runPlan(plan, MPI_COMM_SELF);
  };
  const auto content = [&]() {
    const std::ifstream stream(path);
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
  };

  write("a longer text of an earlier run", false);
  write("short", false);
  CHECK(content() == "short");

  write(", continued", true);
  CHECK(content() == "short, continued");
}

} // namespace seissol::unit_test::runfiles
