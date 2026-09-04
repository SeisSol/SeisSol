// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Instance/Mesh/Xdmf.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "WriterHarness.t.h"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io;

std::string readFile(const std::string& path) {
  const std::ifstream stream(path);
  REQUIRE(stream.good());
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
}

//! Fills the four corners of tetrahedron @p index, as the output does for degree 0.
void projectTestCell(double* target, std::size_t index) {
  std::copy_n(&unit_test::io::TestVertices[index][0][0], 4 * 3, target);
}

} // namespace

// ---------------------------------------------------------------------------
// VTKHDF: write a small mesh and read the file back
// ---------------------------------------------------------------------------

TEST_CASE("IO/VtkHdf: the written file is a well-formed unstructured grid" *
          doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  constexpr std::size_t Cells = 2;
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);

  instance::mesh::VtkHdfWriter vtk(
      "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, false, 0);
  vtk.addPointProjector(projectTestCell);
  vtk.addCellData<std::int32_t>(
      "clustering", {}, true, [](std::int32_t* target, std::size_t index) {
        target[0] = static_cast<std::int32_t>(10 + index);
      });
  vtk.addCellData<double>("v1", {}, false, [](double* target, std::size_t index) {
    target[0] = 0.5 * static_cast<double>(index);
  });

  auto plan = vtk.makeWriter()(dir.prefix(), 0, 0.0);
  unit_test::io::runPlan(plan, MPI_COMM_SELF);

  const auto path = dir.prefix() + "-volume-0.vtkhdf";
  REQUIRE(std::filesystem::exists(path));

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(path);
  hdf5.openGroup("VTKHDF");

  const auto version = hdf5.readAttribute<std::int64_t>("Version");
  REQUIRE(version.size() == 2);
  CHECK(version[0] == 1);

  CHECK(hdf5.readData<std::int64_t>("NumberOfCells").at(0) == static_cast<std::int64_t>(Cells));
  CHECK(hdf5.readData<std::int64_t>("NumberOfPoints").at(0) ==
        static_cast<std::int64_t>(Cells * pointsPerCell));
  CHECK(hdf5.readData<std::int64_t>("NumberOfConnectivityIds").at(0) ==
        static_cast<std::int64_t>(Cells * pointsPerCell));

  // one offset per cell plus the closing one; the last must equal the connectivity length
  const auto offsets = hdf5.readData<std::int64_t>("Offsets");
  REQUIRE(offsets.size() == Cells + 1);
  for (std::size_t i = 0; i < offsets.size(); ++i) {
    CHECK(offsets[i] == static_cast<std::int64_t>(i * pointsPerCell));
  }

  const auto connectivity = hdf5.readData<std::int64_t>("Connectivity");
  REQUIRE(connectivity.size() == Cells * pointsPerCell);
  for (std::size_t i = 0; i < connectivity.size(); ++i) {
    CHECK(connectivity[i] == static_cast<std::int64_t>(i));
  }

  const auto types = hdf5.readData<std::uint8_t>("Types");
  REQUIRE(types.size() == Cells);
  for (const auto type : types) {
    CHECK(type == instance::geometry::vtkType(instance::geometry::Shape::Tetrahedron));
  }

  // the points come back in the order the projector wrote them
  const auto points = hdf5.readData<double>("Points");
  REQUIRE(points.size() == Cells * pointsPerCell * 3);
  for (std::size_t cell = 0; cell < Cells; ++cell) {
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      for (std::size_t d = 0; d < 3; ++d) {
        CHECK(points[(cell * pointsPerCell + vertex) * 3 + d] ==
              doctest::Approx(unit_test::io::TestVertices[cell][vertex][d]));
      }
    }
  }

  hdf5.openGroup("CellData");
  const auto clustering = hdf5.readData<std::int32_t>("clustering");
  REQUIRE(clustering.size() == Cells);
  CHECK(clustering[0] == 10);
  CHECK(clustering[1] == 11);

  const auto v1 = hdf5.readData<double>("v1");
  REQUIRE(v1.size() == Cells);
  CHECK(v1[0] == doctest::Approx(0.0));
  CHECK(v1[1] == doctest::Approx(0.5));
  hdf5.closeGroup();

  hdf5.openGroup("FieldData");
  CHECK(hdf5.readData<double>("Time").at(0) == doctest::Approx(0.0));
  hdf5.closeGroup();

  hdf5.closeGroup();
  hdf5.closeFile();
}

TEST_CASE("IO/VtkHdf: every scheduled write is a self-contained file" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  constexpr std::size_t Cells = 2;

  instance::mesh::VtkHdfWriter vtk(
      "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, true, 0);
  vtk.addPointProjector(projectTestCell);
  double value = 0;
  vtk.addCellData<double>("v1", {}, false, [&value](double* target, std::size_t index) {
    target[0] = value + static_cast<double>(index);
  });

  auto plan = vtk.makeWriter();
  const std::vector<double> times{0.0, 0.25};
  for (std::size_t step = 0; step < times.size(); ++step) {
    value = times[step];
    auto write = plan(dir.prefix(), step, times[step]);
    unit_test::io::runPlan(write, MPI_COMM_SELF);
  }

  // NOTE: the file name carries the step counter even for a time series, so each step produces
  // its own file rather than appending into one. Once WriterGroup::Monolith writes a single file,
  // this is the place to check that Values and CellData grow by one step per write instead.
  for (std::size_t step = 0; step < times.size(); ++step) {
    const auto path = dir.prefix() + "-volume-" + std::to_string(step) + ".vtkhdf";
    REQUIRE_MESSAGE(std::filesystem::exists(path), "missing ", path);

    reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
    hdf5.openFile(path);
    hdf5.openGroup("VTKHDF");

    // a time series announces VTKHDF 2.0 and carries the step bookkeeping
    const auto version = hdf5.readAttribute<std::int64_t>("Version");
    REQUIRE(version.size() == 2);
    CHECK(version[0] == 2);
    CHECK(hdf5.readData<double>("Values").at(0) == doctest::Approx(times[step]));
    CHECK(hdf5.readData<std::uint64_t>("PartOffsets").size() == 1);

    // the geometry is complete in every file
    CHECK(hdf5.readData<std::int64_t>("Connectivity").size() == Cells * 4);

    hdf5.openGroup("CellData");
    const auto v1 = hdf5.readData<double>("v1");
    REQUIRE(v1.size() == Cells);
    CHECK(v1[0] == doctest::Approx(times[step]));
    CHECK(v1[1] == doctest::Approx(times[step] + 1.0));
    hdf5.closeGroup();

    hdf5.closeGroup();
    hdf5.closeFile();
  }
}

TEST_CASE("IO/VtkHdf: a compressed file holds the same data" * doctest::test_suite("io")) {
  // Enough cells that the dataset does not fit into a single chunk, so that the chunking is
  // exercised rather than bypassed.
  constexpr std::size_t Cells = 200000;
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);

  const auto write = [&](const unit_test::io::TempDir& dir, std::int32_t compress) {
    instance::mesh::VtkHdfWriter vtk(
        "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, false, compress);
    vtk.addPointProjector([](double* target, std::size_t index) {
      for (std::size_t vertex = 0; vertex < 4; ++vertex) {
        target[vertex * 3 + 0] = static_cast<double>(index);
        target[vertex * 3 + 1] = static_cast<double>(vertex);
        target[vertex * 3 + 2] = 0.0;
      }
    });
    vtk.addCellData<double>("v1", {}, false, [](double* target, std::size_t index) {
      target[0] = static_cast<double>(index % 3);
    });
    auto plan = vtk.makeWriter()(dir.prefix(), 0, 0.0);
    unit_test::io::runPlan(plan, MPI_COMM_SELF);
    return dir.prefix() + "-volume-0.vtkhdf";
  };

  const unit_test::io::TempDir plainDir;
  const unit_test::io::TempDir compressedDir;
  const auto plain = write(plainDir, 0);
  const auto compressed = write(compressedDir, 6);

  CHECK(std::filesystem::file_size(compressed) < std::filesystem::file_size(plain));

  // compression must not change what is stored
  for (const auto& path : {plain, compressed}) {
    reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
    hdf5.openFile(path);
    hdf5.openGroup("VTKHDF");

    const auto connectivity = hdf5.readData<std::int64_t>("Connectivity");
    REQUIRE(connectivity.size() == Cells * pointsPerCell);
    CHECK(connectivity.front() == 0);
    CHECK(connectivity.back() == static_cast<std::int64_t>(Cells * pointsPerCell - 1));

    hdf5.openGroup("CellData");
    const auto values = hdf5.readData<double>("v1");
    REQUIRE(values.size() == Cells);
    for (std::size_t index = 0; index < Cells; index += 4096) {
      CHECK(values[index] == doctest::Approx(static_cast<double>(index % 3)));
    }
    hdf5.closeGroup();

    hdf5.closeGroup();
    hdf5.closeFile();
  }
}

TEST_CASE("IO/VtkHdf: an incremental snapshot links to the constant data" *
          doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  constexpr std::size_t Cells = 2;
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);

  instance::mesh::VtkHdfWriter vtk(
      "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, false, 0, true);
  vtk.addPointProjector(projectTestCell);
  vtk.addCellData<std::int32_t>(
      "clustering", {}, true, [](std::int32_t* target, std::size_t index) {
        target[0] = static_cast<std::int32_t>(10 + index);
      });
  double value = 0;
  vtk.addCellData<double>("v1", {}, false, [&value](double* target, std::size_t index) {
    target[0] = value + static_cast<double>(index);
  });

  auto plan = vtk.makeWriter();
  const std::vector<double> times{0.0, 0.25};
  for (std::size_t step = 0; step < times.size(); ++step) {
    value = times[step];
    auto write = plan(dir.prefix(), step, times[step]);
    unit_test::io::runPlan(write, MPI_COMM_SELF);
  }

  // the geometry and the constant cell data go into one file of their own
  const auto constFile = dir.prefix() + "-volume-const.vtkhdf";
  REQUIRE(std::filesystem::exists(constFile));
  const auto constSize = std::filesystem::file_size(constFile);

  for (std::size_t step = 0; step < times.size(); ++step) {
    const auto path = dir.prefix() + "-volume-" + std::to_string(step) + ".vtkhdf";
    REQUIRE_MESSAGE(std::filesystem::exists(path), "missing ", path);

    // a snapshot only holds what changes, so it stays smaller than the constant part
    CHECK(std::filesystem::file_size(path) < constSize);

    // ... and the links make it look complete to a reader
    reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
    hdf5.openFile(path);
    hdf5.openGroup("VTKHDF");

    const auto connectivity = hdf5.readData<std::int64_t>("Connectivity");
    REQUIRE(connectivity.size() == Cells * pointsPerCell);
    CHECK(connectivity.back() == static_cast<std::int64_t>(Cells * pointsPerCell - 1));

    const auto points = hdf5.readData<double>("Points");
    REQUIRE(points.size() == Cells * pointsPerCell * 3);
    CHECK(points[0] == doctest::Approx(unit_test::io::TestVertices[0][0][0]));

    hdf5.openGroup("CellData");
    const auto clustering = hdf5.readData<std::int32_t>("clustering");
    REQUIRE(clustering.size() == Cells);
    CHECK(clustering[0] == 10);

    // the changing data is in the snapshot itself, and differs per step
    const auto v1 = hdf5.readData<double>("v1");
    REQUIRE(v1.size() == Cells);
    CHECK(v1[0] == doctest::Approx(times[step]));
    CHECK(v1[1] == doctest::Approx(times[step] + 1.0));
    hdf5.closeGroup();

    hdf5.closeGroup();
    hdf5.closeFile();
  }
}

// ---------------------------------------------------------------------------
// Xdmf: the declared dimensions have to match the payload that is written
// ---------------------------------------------------------------------------

TEST_CASE("IO/Xdmf: the XML describes the payload that was written" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  constexpr std::size_t Cells = 2;
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);

  instance::mesh::XdmfWriter xdmf(
      "volume", Cells, instance::geometry::Shape::Tetrahedron, 0, false, 0);
  xdmf.addPointProjector(projectTestCell);
  xdmf.addCellData<double>("v1", {}, false, [](double* target, std::size_t index) {
    target[0] = 0.5 * static_cast<double>(index);
  });

  auto plan = xdmf.makeWriter()(dir.prefix(), 0, 0.0);
  unit_test::io::runPlan(plan, MPI_COMM_SELF);

  const auto xml = readFile(dir.prefix() + "-volume.xdmf");

  // the topology and geometry sizes are the first thing a reader trusts
  CHECK(xml.find("TopologyType=\"Tetrahedron\"") != std::string::npos);
  CHECK(xml.find("NumberOfElements=\"" + std::to_string(Cells) + "\"") != std::string::npos);
  CHECK(xml.find("Dimensions=\"" + std::to_string(Cells * pointsPerCell) + " 3\"") !=
        std::string::npos);
  CHECK(xml.find("Dimensions=\"" + std::to_string(Cells) + " " + std::to_string(pointsPerCell) +
                 "\"") != std::string::npos);
  CHECK(xml.find("Name=\"v1\"") != std::string::npos);
  CHECK(xml.find("Center=\"Cell\"") != std::string::npos);

  // the payload the XML points at has to exist and hold what was announced
  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(dir.prefix() + "-volume-data.h5");
  // dataset0 is the topology, dataset1 the geometry: one row per point, three coordinates each
  CHECK(hdf5.dataCount("dataset0") == Cells);
  CHECK(hdf5.dataRowSize("dataset0") == pointsPerCell);
  CHECK(hdf5.dataCount("dataset1") == Cells * pointsPerCell);
  CHECK(hdf5.dataRowSize("dataset1") == 3);
  const auto geometry = hdf5.readData<double>("dataset1");
  REQUIRE(geometry.size() == Cells * pointsPerCell * 3);
  for (std::size_t cell = 0; cell < Cells; ++cell) {
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      for (std::size_t d = 0; d < 3; ++d) {
        CHECK(geometry[(cell * pointsPerCell + vertex) * 3 + d] ==
              doctest::Approx(unit_test::io::TestVertices[cell][vertex][d]));
      }
    }
  }
  hdf5.closeFile();
}

} // namespace seissol::unit_test
