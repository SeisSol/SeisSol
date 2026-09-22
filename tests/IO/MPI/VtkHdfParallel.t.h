// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Mesh/VtkHdf.h"
#include "IO/Instance/Point/Grouping.h"
#include "IO/Instance/Point/Hdf5Table.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/WriterHarness.t.h"
#include "Parallel/MPI.h"

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

/**
 * With shared points the number a rank contributes is no longer its cell count times the corners
 * of a cell, so the point offsets come out of their own scan rather than out of the element one.
 * That scan is what this checks.
 */
TEST_CASE("IO/VtkHdf: shared points are offset per rank" * doctest::test_suite("io")) {
  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  REQUIRE(size > 1);

  // rank r writes r cells, all four corners of each at the same place, so it contributes exactly
  // one point -- and none at all when it has no cells
  const auto localCells = rank;
  const auto globalCells = (size * (size - 1)) / 2;
  const auto contributingRanks = size - 1;

  const SharedTempDir dir;
  instance::geometry::WriterConfig config;
  config.order = 0;
  config.format = instance::geometry::WriterFormat::Vtk;

  instance::geometry::GeometryWriter geometry(
      "volume",
      localCells,
      instance::geometry::Shape::Tetrahedron,
      config,
      1,
      [=](double* target, std::size_t /*cell*/, std::size_t /*subcell*/) {
        for (std::size_t corner = 0; corner < 4; ++corner) {
          target[corner * 3 + 0] = static_cast<double>(rank);
          target[corner * 3 + 1] = 0.0;
          target[corner * 3 + 2] = 0.0;
        }
      });
  geometry.addGeometryOutput<double>(
      "v1", {}, false, [=](double* target, std::size_t cell, std::size_t /*subcell*/) {
        target[0] = static_cast<double>(rank * 1000 + cell);
      });

  auto plan = geometry.makeWriter()(dir.prefix(), 0, 0.0);
  unit_test::io::runPlan(plan, Mpi::mpi.comm());
  MPI_Barrier(Mpi::mpi.comm());

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(dir.prefix() + "-volume-0.vtkhdf");
  hdf5.openGroup("VTKHDF");

  CHECK(hdf5.readData<std::int64_t>("NumberOfCells").at(0) ==
        static_cast<std::int64_t>(globalCells));
  // one point per contributing rank, not four per cell
  CHECK(hdf5.readData<std::int64_t>("NumberOfPoints").at(0) ==
        static_cast<std::int64_t>(contributingRanks));
  CHECK(hdf5.readData<std::int64_t>("NumberOfConnectivityIds").at(0) ==
        static_cast<std::int64_t>(globalCells * 4));

  // every cell points at the one point of the rank that wrote it
  const auto connectivity = hdf5.readData<std::int64_t>("Connectivity");
  const auto points = hdf5.readData<double>("Points");
  REQUIRE(connectivity.size() == globalCells * 4);
  REQUIRE(points.size() == contributingRanks * 3);
  for (std::size_t entry = 0; entry < connectivity.size(); ++entry) {
    const auto point = connectivity[entry];
    REQUIRE(point >= 0);
    REQUIRE(static_cast<std::size_t>(point) < contributingRanks);
    // rank 0 contributes nothing, so the point of rank r is the (r-1)th
    CHECK(points[point * 3] == doctest::Approx(static_cast<double>(point + 1)));
  }

  hdf5.closeGroup();
  hdf5.closeFile();
}

/**
 * The groups and their order have to be the same on every rank, including on ranks that hold no
 * point of a group at all -- they still take part in declaring its dataset. And the renumbering
 * has to hand every rank one contiguous range per group, so that the group stays expressible as a
 * single distributed dimension.
 */
TEST_CASE("IO/Grouping: the groups agree across ranks" * doctest::test_suite("io")) {
  using namespace seissol::io::instance::point;
  using seissol::io::datatype::inferDatatype;

  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  if (size < 3) {
    // with fewer ranks the setup below cannot produce both quantity sets at once
    return;
  }

  const std::vector<TableQuantity> elastic{{"v1", inferDatatype<double>()},
                                           {"v2", inferDatatype<double>()}};
  const std::vector<TableQuantity> poroelastic{{"v1", inferDatatype<double>()},
                                               {"v2", inferDatatype<double>()},
                                               {"p", inferDatatype<double>()}};

  // only the last rank holds poroelastic points, and rank 0 holds none at all
  std::vector<std::vector<TableQuantity>> points;
  for (std::size_t point = 0; point < rank; ++point) {
    points.push_back(rank + 1 == size ? poroelastic : elastic);
  }

  const auto grouping = groupPoints(points, Mpi::mpi.comm());

  // both groups exist everywhere, in the same order
  REQUIRE(grouping.groupCount() == 2);
  const auto poroGroup = grouping.quantities[0].size() == 3 ? 0 : 1;
  const auto elasticGroup = 1 - poroGroup;
  CHECK(grouping.quantities[elasticGroup].size() == 2);

  const auto elasticTotal = ((size - 1) * (size - 2)) / 2;
  CHECK(grouping.globalCount[elasticGroup] == elasticTotal);
  CHECK(grouping.globalCount[poroGroup] == size - 1);

  // this rank's indices form one contiguous range inside its group
  if (!points.empty()) {
    const auto group = grouping.group.front();
    for (const auto other : grouping.group) {
      CHECK(other == group);
    }
    for (std::size_t point = 1; point < grouping.index.size(); ++point) {
      CHECK(grouping.index[point] == grouping.index[point - 1] + 1);
    }
    CHECK(grouping.index.back() < grouping.globalCount[group]);
  }

  // ... and the ranges of the ranks cover the group exactly once
  std::vector<int> covered(grouping.globalCount[elasticGroup] + grouping.globalCount[poroGroup], 0);
  std::vector<int> mine(covered.size(), 0);
  for (std::size_t point = 0; point < points.size(); ++point) {
    const auto flat = grouping.index[point] +
                      (grouping.group[point] == poroGroup ? grouping.globalCount[elasticGroup] : 0);
    ++mine[flat];
  }
  MPI_Allreduce(mine.data(),
                covered.data(),
                static_cast<int>(covered.size()),
                MPI_INT,
                MPI_SUM,
                Mpi::mpi.comm());
  for (const auto count : covered) {
    CHECK(count == 1);
  }
}

TEST_CASE("IO/Hdf5Table: a rank without points of a table still writes with it" *
          doctest::test_suite("io")) {
  using namespace seissol::io::instance::point;

  const SharedTempDir dir;
  const auto rank = Mpi::mpi.rank();
  const auto size = Mpi::mpi.size();

  const auto quantity = [](const std::string& name) {
    return TableQuantity{name, datatype::inferDatatype<double>()};
  };
  const std::vector<TableQuantity> small{quantity("v1"), quantity("v2")};
  const std::vector<TableQuantity> large{quantity("v1"), quantity("v2"), quantity("p")};

  // only rank 0 holds a point of the wider set, so every other rank takes part in declaring and
  // appending to a table it has nothing to put in
  std::vector<std::vector<TableQuantity>> pointQuantities{small};
  if (rank == 0) {
    pointQuantities.push_back(large);
  }

  Hdf5Table table("points", pointQuantities, Mpi::mpi.comm(), 4);
  const auto& grouping = table.grouping();
  REQUIRE(grouping.groupCount() == 2);

  const auto smallGroup = grouping.group[0];
  const auto largeGroup = smallGroup == 0 ? 1 : 0;
  CHECK(table.localPointCount(smallGroup) == 1);
  CHECK(table.localPointCount(largeGroup) == (rank == 0 ? 1U : 0U));
  CHECK(grouping.globalCount[smallGroup] == static_cast<std::size_t>(size));
  CHECK(grouping.globalCount[largeGroup] == 1);

  auto plan = table.makeWriter();
  const std::vector<std::size_t> samplesPerWrite{2, 3};
  std::size_t written = 0;
  for (std::size_t step = 0; step < samplesPerWrite.size(); ++step) {
    const auto samples = samplesPerWrite[step];
    for (std::size_t group = 0; group < grouping.groupCount(); ++group) {
      auto* storage = reinterpret_cast<double*>(table.prepare(group, samples));
      const auto points = table.localPointCount(group);
      const auto components = table.sampleSize(group) / sizeof(double);
      for (std::size_t sample = 0; sample < samples; ++sample) {
        for (std::size_t point = 0; point < points; ++point) {
          for (std::size_t component = 0; component < components; ++component) {
            storage[(sample * points + point) * components + component] =
                1000.0 * rank + 10.0 * static_cast<double>(written + sample) +
                static_cast<double>(component);
          }
        }
      }
    }
    auto write = plan(dir.prefix(), step, static_cast<double>(step));
    unit_test::io::runPlan(write, Mpi::mpi.comm());
    written += samples;
  }

  MPI_Barrier(Mpi::mpi.comm());

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(dir.prefix() + "-points.h5");
  hdf5.openGroup("points");

  // the narrow table holds a row per rank and sample, the wide one only the row of rank 0
  const auto narrow = hdf5.readData<double>("group" + std::to_string(smallGroup));
  CHECK(narrow.size() == written * static_cast<std::size_t>(size) * 2);
  const auto wide = hdf5.readData<double>("group" + std::to_string(largeGroup));
  CHECK(wide.size() == written * 3);

  // and the row of a rank sits where the grouping said it would
  for (std::size_t sample = 0; sample < written; ++sample) {
    for (int other = 0; other < size; ++other) {
      const auto base = (sample * static_cast<std::size_t>(size) + other) * 2;
      CHECK(narrow[base] == doctest::Approx(1000.0 * other + 10.0 * static_cast<double>(sample)));
    }
  }

  hdf5.closeGroup();
  hdf5.closeFile();
}

} // namespace seissol::unit_test
