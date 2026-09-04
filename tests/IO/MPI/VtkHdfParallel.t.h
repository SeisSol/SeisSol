// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "Parallel/MPI.h"
#include "WriterHarness.t.h"

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io;

/**
 * A temporary directory that every rank agrees on: rank 0 creates it and broadcasts the name.
 */
struct SharedTempDir {
  std::string path;

  SharedTempDir() {
    std::vector<char> buffer(256, '\0');
    if (Mpi::mpi.rank() == 0) {
      std::string name = "/tmp/seissoliompitestXXXXXX";
      const char* created = mkdtemp(name.data());
      REQUIRE(created != nullptr);
      std::copy(name.begin(), name.end(), buffer.begin());
    }
    MPI_Bcast(buffer.data(), static_cast<int>(buffer.size()), MPI_CHAR, 0, Mpi::mpi.comm());
    path = std::string(buffer.data());
    // the ranks may not share a file system view instantly; make sure the directory is there
    // everywhere before anyone opens a file in it collectively
    std::error_code error;
    std::filesystem::create_directories(path, error);
    MPI_Barrier(Mpi::mpi.comm());
  }

  ~SharedTempDir() {
    MPI_Barrier(Mpi::mpi.comm());
    if (Mpi::mpi.rank() == 0) {
      std::error_code error;
      std::filesystem::remove_all(path, error);
    }
  }

  SharedTempDir(const SharedTempDir&) = delete;
  SharedTempDir(SharedTempDir&&) = delete;
  auto operator=(const SharedTempDir&) -> SharedTempDir& = delete;
  auto operator=(SharedTempDir&&) -> SharedTempDir& = delete;

  [[nodiscard]] std::string prefix() const { return path + "/out"; }
};

} // namespace

/**
 * The cell and point offsets of the VTKHDF output come out of an MPI_Exscan over the local element
 * counts, and the closing entry of Offsets is only contributed by the last rank. Both break in
 * ways a single-rank test cannot see, in particular when a partition is empty -- which happens in
 * practice as soon as an output region filter excludes a whole rank.
 */
TEST_CASE("IO/VtkHdf: unequal partitions produce one consistent grid" * doctest::test_suite("io")) {
  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  REQUIRE(size > 1);

  // rank r contributes r cells, so rank 0 stays empty
  const auto localCells = rank;
  const auto globalCells = (size * (size - 1)) / 2;
  const auto offset = (rank * (rank - 1)) / 2;
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);

  const SharedTempDir dir;

  instance::mesh::VtkHdfWriter vtk(
      "volume", localCells, instance::geometry::Shape::Tetrahedron, 0, false, 0);
  vtk.addPointProjector([=](double* target, std::size_t index) {
    // a degenerate but unique cell per global index, so that the order can be checked
    const auto global = static_cast<double>(offset + index);
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      target[vertex * 3 + 0] = global;
      target[vertex * 3 + 1] = static_cast<double>(vertex);
      target[vertex * 3 + 2] = 0.0;
    }
  });
  // the value encodes which rank wrote the cell and where it sat locally
  vtk.addCellData<double>("v1", {}, false, [=](double* target, std::size_t index) {
    target[0] = 1000.0 * static_cast<double>(rank) + static_cast<double>(index);
  });

  auto plan = vtk.makeWriter()(dir.prefix(), 0, 0.0);
  unit_test::io::runPlan(plan, Mpi::mpi.comm());

  MPI_Barrier(Mpi::mpi.comm());

  const auto path = dir.prefix() + "-volume-0.vtkhdf";
  REQUIRE(std::filesystem::exists(path));

  // read the whole file on every rank, so the assertions below see the global picture
  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(path);
  hdf5.openGroup("VTKHDF");

  CHECK(hdf5.readData<std::int64_t>("NumberOfCells").at(0) ==
        static_cast<std::int64_t>(globalCells));
  CHECK(hdf5.readData<std::int64_t>("NumberOfPoints").at(0) ==
        static_cast<std::int64_t>(globalCells * pointsPerCell));

  // exactly one closing entry overall, not one per rank
  const auto offsets = hdf5.readData<std::int64_t>("Offsets");
  REQUIRE(offsets.size() == globalCells + 1);
  for (std::size_t i = 0; i < offsets.size(); ++i) {
    CHECK(offsets[i] == static_cast<std::int64_t>(i * pointsPerCell));
  }

  // the connectivity is the concatenation of the per-rank blocks, without gaps or overlaps
  const auto connectivity = hdf5.readData<std::int64_t>("Connectivity");
  REQUIRE(connectivity.size() == globalCells * pointsPerCell);
  for (std::size_t i = 0; i < connectivity.size(); ++i) {
    CHECK(connectivity[i] == static_cast<std::int64_t>(i));
  }

  CHECK(hdf5.readData<std::uint8_t>("Types").size() == globalCells);

  // every cell sits where its rank's offset says it should
  const auto points = hdf5.readData<double>("Points");
  REQUIRE(points.size() == globalCells * pointsPerCell * 3);

  hdf5.openGroup("CellData");
  const auto values = hdf5.readData<double>("v1");
  REQUIRE(values.size() == globalCells);
  for (std::size_t other = 0; other < size; ++other) {
    const auto otherOffset = (other * (other - 1)) / 2;
    for (std::size_t index = 0; index < other; ++index) {
      const auto expected = 1000.0 * static_cast<double>(other) + static_cast<double>(index);
      CHECK(values[otherOffset + index] == doctest::Approx(expected));
      CHECK(points[(otherOffset + index) * pointsPerCell * 3] ==
            doctest::Approx(static_cast<double>(otherOffset + index)));
    }
  }
  hdf5.closeGroup();

  hdf5.closeGroup();
  hdf5.closeFile();
}

/**
 * A dataset with a filter is, in parallel HDF5, gathered chunk by chunk onto a single rank before
 * it is compressed. That path is separate from the plain one and only runs on more than one rank.
 */
TEST_CASE("IO/VtkHdf: compression survives unequal partitions" * doctest::test_suite("io")) {
  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  REQUIRE(size > 1);

  // enough cells to span several chunks, and again an empty rank 0
  const auto localCells = rank * 100000;
  const auto globalCells = ((size * (size - 1)) / 2) * 100000;
  const auto offset = ((rank * (rank - 1)) / 2) * 100000;
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);

  const SharedTempDir dir;

  instance::mesh::VtkHdfWriter vtk(
      "volume", localCells, instance::geometry::Shape::Tetrahedron, 0, false, 6);
  vtk.addPointProjector([=](double* target, std::size_t index) {
    const auto global = static_cast<double>(offset + index);
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      target[vertex * 3 + 0] = global;
      target[vertex * 3 + 1] = static_cast<double>(vertex);
      target[vertex * 3 + 2] = 0.0;
    }
  });
  vtk.addCellData<double>("v1", {}, false, [=](double* target, std::size_t index) {
    target[0] = static_cast<double>((offset + index) % 7);
  });

  auto plan = vtk.makeWriter()(dir.prefix(), 0, 0.0);
  unit_test::io::runPlan(plan, Mpi::mpi.comm());

  MPI_Barrier(Mpi::mpi.comm());

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(dir.prefix() + "-volume-0.vtkhdf");
  hdf5.openGroup("VTKHDF");

  const auto connectivity = hdf5.readData<std::int64_t>("Connectivity");
  REQUIRE(connectivity.size() == globalCells * pointsPerCell);
  CHECK(connectivity.back() == static_cast<std::int64_t>(globalCells * pointsPerCell - 1));

  hdf5.openGroup("CellData");
  const auto values = hdf5.readData<double>("v1");
  REQUIRE(values.size() == globalCells);
  for (std::size_t index = 0; index < globalCells; index += 4096) {
    CHECK(values[index] == doctest::Approx(static_cast<double>(index % 7)));
  }
  hdf5.closeGroup();

  hdf5.closeGroup();
  hdf5.closeFile();
}

} // namespace seissol::unit_test
