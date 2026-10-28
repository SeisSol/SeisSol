// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Instance/Mesh/Xdmf.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/Writer/File/Hdf5Writer.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Dimension.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "WriterHarness.t.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
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

//! The backups an earlier run left for the files whose names start with @p stem .
inline std::vector<std::string> backupsOf(const std::string& directory, const std::string& stem) {
  std::vector<std::string> found;
  for (const auto& entry : std::filesystem::directory_iterator(directory)) {
    const auto name = entry.path().filename().string();
    if (name.rfind(stem, 0) == 0 && name.find(".bak_") != std::string::npos) {
      found.push_back(entry.path().string());
    }
  }
  return found;
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

  // a run resuming from a checkpoint writes a file of its own, and keeps the earlier one
  writer::file::RunFiles resumed;
  resumed.resumed = true;
  auto next = appendPlan(path, 4.0);
  unit_test::io::runPlan(next, MPI_COMM_SELF, &resumed);
  CHECK(readBack(path, "series") == std::vector<double>{4.0});
  const auto backups = backupsOf(dir.path, "out-series");
  REQUIRE(backups.size() == 1);
  CHECK(readBack(backups[0], "series") == std::vector<double>{1.0, 2.0, 3.0});

  // and a new run starts it over, without a backup
  writer::file::RunFiles fresh;
  auto first = appendPlan(path, 5.0);
  unit_test::io::runPlan(first, MPI_COMM_SELF, &fresh);
  CHECK(readBack(path, "series") == std::vector<double>{5.0});
  CHECK(backupsOf(dir.path, "out-series").size() == 1);
}

TEST_CASE("IO/RunFiles: a resumed time series is a file of its own" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  constexpr std::size_t Cells = 2;
  const std::vector<double> times{0.0, 0.25, 0.5, 0.75, 1.0};

  // the run before the checkpoint writes steps 0 to 2, the resumed one continues at step 3
  const auto run = [&](std::size_t from, std::size_t to, writer::file::RunFiles& files) {
    instance::mesh::VtkHdfWriter vtk(
        "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, true, 0);
    vtk.addPointProjector([](double* target, std::size_t index) {
      std::copy_n(&unit_test::io::TestVertices[index][0][0], 4 * 3, target);
    });
    double value = 0;
    vtk.addCellData<double>("v1", {}, false, [&value](double* target, std::size_t index) {
      target[0] = value + static_cast<double>(index);
    });
    auto plan = vtk.makeWriter();
    for (std::size_t counter = from; counter < to; ++counter) {
      value = times[counter];
      auto write = plan(dir.prefix(), counter, times[counter]);
      unit_test::io::runPlan(write, MPI_COMM_SELF, &files);
    }
  };
  writer::file::RunFiles before;
  run(0, 3, before);
  writer::file::RunFiles resumed;
  resumed.resumed = true;
  run(3, 5, resumed);

  REQUIRE(backupsOf(dir.path, "out-volume").size() == 1);

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(dir.prefix() + "-volume.vtkhdf");
  hdf5.openGroup("VTKHDF");
  // the file carries the mesh, and the two steps of the resumed run, counted from its start
  CHECK(hdf5.readData<std::int64_t>("NumberOfCells").size() == 1);
  CHECK(hdf5.readData<std::int64_t>("Connectivity").size() == Cells * 4);
  hdf5.openGroup("CellData");
  const auto v1 = hdf5.readData<double>("v1");
  REQUIRE(v1.size() == 2 * Cells);
  CHECK(v1[0] == doctest::Approx(times[3]));
  CHECK(v1[Cells] == doctest::Approx(times[4]));
  hdf5.closeGroup();
  hdf5.openGroup("Steps");
  CHECK(hdf5.readAttributeScalar<std::int64_t>("NSteps") == 2);
  CHECK(hdf5.readData<double>("Values") == std::vector<double>{times[3], times[4]});
  hdf5.openGroup("CellDataOffsets");
  CHECK(hdf5.readData<std::uint64_t>("v1") == std::vector<std::uint64_t>{0, Cells});
  hdf5.closeGroup();
  hdf5.closeGroup();
  // while the index of an output keeps counting across the runs
  hdf5.openGroup("FieldData");
  CHECK(hdf5.readData<std::size_t>("Index") == std::vector<std::size_t>{3, 4});
  hdf5.closeGroup();
  hdf5.closeGroup();
  hdf5.closeFile();
}

TEST_CASE("IO/RunFiles: a resumed Xdmf output describes its own files" *
          doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  constexpr std::size_t Cells = 2;

  const auto run = [&](std::size_t from, std::size_t to, writer::file::RunFiles& files) {
    // the binary backend, whose payload files are plain byte streams
    instance::mesh::XdmfWriter xdmf(
        "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, true, 0);
    xdmf.addPointProjector([](double* target, std::size_t index) {
      std::copy_n(&unit_test::io::TestVertices[index][0][0], 4 * 3, target);
    });
    std::size_t current = 0;
    xdmf.addCellData<double>("v1", {}, false, [&current](double* target, std::size_t index) {
      target[0] = 10.0 * static_cast<double>(current) + static_cast<double>(index);
    });
    auto plan = xdmf.makeWriter();
    for (std::size_t counter = from; counter < to; ++counter) {
      current = counter;
      auto write = plan(dir.prefix(), counter, static_cast<double>(counter));
      unit_test::io::runPlan(write, MPI_COMM_SELF, &files);
    }
  };
  writer::file::RunFiles before;
  run(0, 3, before);
  writer::file::RunFiles resumed;
  resumed.resumed = true;
  run(3, 5, resumed);

  // the mesh is described again, and the two steps of the resumed run lie at its start
  const auto xml = [&]() {
    const std::ifstream stream(dir.prefix() + "-volume.xdmf");
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
  }();
  CHECK(xml.find("TopologyType=\"Tetrahedron\"") != std::string::npos);
  CHECK(xml.find("GeometryType=\"XYZ\"") != std::string::npos);
  CHECK(xml.find("step-3") != std::string::npos);
  CHECK(xml.find("step-4") != std::string::npos);
  CHECK(xml.find("step-2") == std::string::npos);

  // every file of the output was kept from the run before
  CHECK(backupsOf(dir.path, "out-volume").size() >= 2);
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
