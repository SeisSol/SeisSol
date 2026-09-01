// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Geometry/Geometry.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Instructions/Instruction.h"
#include "IO/Writer/Writer.h"
#include "WriterHarness.t.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io;
using namespace seissol::io::writer;

//! The in-file path an instruction writes to, plus its dataset name.
std::string target(const std::shared_ptr<instructions::WriteInstruction>& instruction) {
  if (auto* data = dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get()); data != nullptr) {
    return data->location.infilePath() + "/" + data->name;
  }
  if (auto* attribute = dynamic_cast<instructions::Hdf5AttributeWrite*>(instruction.get());
      attribute != nullptr) {
    return attribute->location.infilePath() + "@" + attribute->name;
  }
  return {};
}

std::vector<std::string> targets(writer::Writer& writer) {
  std::vector<std::string> result;
  for (const auto& instruction : writer.getInstructions()) {
    result.push_back(target(instruction));
  }
  return result;
}

bool contains(const std::vector<std::string>& haystack, const std::string& needle) {
  return std::find(haystack.begin(), haystack.end(), needle) != haystack.end();
}

//! The value of an inline int64 dataset, so that the announced global counts can be checked.
std::int64_t inlineValue(writer::Writer& writer, const std::string& path) {
  for (const auto& instruction : writer.getInstructions()) {
    auto* data = dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get());
    if (data == nullptr || target(instruction) != path) {
      continue;
    }
    const auto* source = dynamic_cast<WriteInline*>(data->dataSource.get());
    REQUIRE(source != nullptr);
    std::int64_t value = 0;
    std::copy_n(static_cast<const char*>(source->getLocalPointer()),
                sizeof(value),
                reinterpret_cast<char*>(&value));
    return value;
  }
  FAIL("no inline dataset at ", path);
  return 0;
}

} // namespace

// ---------------------------------------------------------------------------
// Instruction and data source serialisation
// ---------------------------------------------------------------------------

TEST_CASE("IO/Instruction: Hdf5Location round-trips" * doctest::test_suite("io")) {
  const instructions::Hdf5Location location("file.h5", {"VTKHDF", "CellData"}, "partition");

  CHECK(location.file() == "file.h5");
  REQUIRE(location.groups().size() == 2);
  CHECK(location.groups()[0] == "VTKHDF");
  CHECK(location.groups()[1] == "CellData");
  REQUIRE(location.dataset().has_value());
  CHECK(location.dataset().value() == "partition");
  // infilePath() includes the dataset when the location names one
  CHECK(location.infilePath() == "/VTKHDF/CellData/partition");

  auto copy = location;
  instructions::Hdf5Location restored(copy.serialize());
  CHECK(restored.file() == location.file());
  CHECK(restored.groups() == location.groups());
  CHECK(restored.infilePath() == location.infilePath());
}

TEST_CASE("IO/Instruction: common location of two paths" * doctest::test_suite("io")) {
  const instructions::Hdf5Location cellData("file.h5", {"VTKHDF", "CellData"});
  const instructions::Hdf5Location pointData("file.h5", {"VTKHDF", "PointData"});
  const instructions::Hdf5Location other("elsewhere.h5", {"VTKHDF"});

  const auto common = cellData.commonLocation(pointData);
  REQUIRE(common.has_value());
  CHECK(common.value().infilePath() == "/VTKHDF");

  CHECK_FALSE(cellData.commonLocation(other).has_value());
}

TEST_CASE("IO/Instruction: inline data survives serialisation" * doctest::test_suite("io")) {
  const std::vector<std::int64_t> values{3, 1, 4, 1, 5};
  const auto source = WriteInline::createArray<std::int64_t>({values.size()}, values);

  const auto restored = DataSource::deserialize(source->serialize());
  REQUIRE(restored->getLocalSize() == values.size() * sizeof(std::int64_t));
  const auto* restoredValues = static_cast<const std::int64_t*>(restored->getLocalPointer());
  for (std::size_t i = 0; i < values.size(); ++i) {
    CHECK(restoredValues[i] == values[i]);
  }
  CHECK_FALSE(restored->distributed());
}

TEST_CASE("IO/Instruction: a plan survives serialisation" * doctest::test_suite("io")) {
  writer::Writer original;
  original.addInstruction(std::make_shared<instructions::Hdf5DataWrite>(
      instructions::Hdf5Location("file.h5", {"VTKHDF"}),
      "NumberOfCells",
      WriteInline::createArray<std::int64_t>({1}, {7}),
      datatype::inferDatatype<std::int64_t>()));
  original.addInstruction(std::make_shared<instructions::Hdf5AttributeWrite>(
      instructions::Hdf5Location("file.h5", {"VTKHDF"}),
      "Version",
      WriteInline::createArray<std::int64_t>({2}, {1, 0})));

  const auto serialized = original.serialize();
  writer::Writer restored(serialized);

  REQUIRE(restored.getInstructions().size() == original.getInstructions().size());
  CHECK(targets(restored) == targets(original));
  CHECK(restored.serialize() == serialized);
  CHECK(inlineValue(restored, "/VTKHDF/NumberOfCells") == 7);
}

// ---------------------------------------------------------------------------
// The plans the mesh writers produce
// ---------------------------------------------------------------------------

TEST_CASE("IO/Instruction: the VTKHDF plan has the mandatory datasets" *
          doctest::test_suite("io")) {
  // cell output (degree 0) of two tetrahedra
  instance::mesh::VtkHdfWriter vtk(
      "volume", 2, instance::geometry::Shape::Tetrahedron, 0, false, 0);
  vtk.addPointProjector([](double* target, std::size_t index) {
    std::copy_n(&unit_test::io::TestVertices[index][0][0], 4 * 3, target);
  });
  vtk.addCellData<std::int32_t>(
      "partition", {}, true, [](std::int32_t* target, std::size_t /*index*/) { target[0] = 0; });

  auto plan = vtk.makeWriter()("prefix", 0, 0.0);
  const auto paths = targets(plan);

  // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format
  for (const auto* mandatory : {"/VTKHDF@Type",
                                "/VTKHDF@Version",
                                "/VTKHDF/NumberOfCells",
                                "/VTKHDF/NumberOfPoints",
                                "/VTKHDF/NumberOfConnectivityIds",
                                "/VTKHDF/Points",
                                "/VTKHDF/Offsets",
                                "/VTKHDF/Types",
                                "/VTKHDF/Connectivity"}) {
    CHECK_MESSAGE(contains(paths, mandatory), "missing ", mandatory);
  }
  CHECK(contains(paths, "/VTKHDF/CellData/partition"));
  CHECK(contains(paths, "/VTKHDF/FieldData/Time"));
  CHECK(contains(paths, "/VTKHDF/FieldData/Index"));

  // one point per vertex of every cell, no vertex sharing
  const auto pointsPerCell =
      instance::geometry::numPoints(1, instance::geometry::Shape::Tetrahedron);
  CHECK(inlineValue(plan, "/VTKHDF/NumberOfCells") == 2);
  CHECK(inlineValue(plan, "/VTKHDF/NumberOfPoints") ==
        2 * static_cast<std::int64_t>(pointsPerCell));
  CHECK(inlineValue(plan, "/VTKHDF/NumberOfConnectivityIds") ==
        2 * static_cast<std::int64_t>(pointsPerCell));
}

TEST_CASE("IO/Instruction: the time series adds the step bookkeeping" * doctest::test_suite("io")) {
  instance::mesh::VtkHdfWriter vtk("volume", 2, instance::geometry::Shape::Tetrahedron, 0, true, 0);
  vtk.addCellData<double>(
      "v1", {}, false, [](double* target, std::size_t index) { target[0] = index; });

  auto plan = vtk.makeWriter()("prefix", 0, 0.5);
  const auto paths = targets(plan);

  CHECK(contains(paths, "/VTKHDF/Values"));
  CHECK(contains(paths, "/VTKHDF/PartOffsets"));
  CHECK(contains(paths, "/VTKHDF/PointOffsets"));
  CHECK(contains(paths, "/VTKHDF/CellOffsets"));
  CHECK(contains(paths, "/VTKHDF/ConnectivityIdOffsets"));
  CHECK(contains(paths, "/VTKHDF/CellDataOffsets/v1"));

  // non-const data is appended, const data is not
  for (const auto& instruction : plan.getInstructions()) {
    auto* data = dynamic_cast<instructions::Hdf5DataWrite*>(instruction.get());
    if (data == nullptr) {
      continue;
    }
    if (target(instruction) == "/VTKHDF/CellData/v1") {
      CHECK(data->append);
    }
    if (target(instruction) == "/VTKHDF/Connectivity") {
      CHECK_FALSE(data->append);
    }
  }
}

TEST_CASE("IO/Instruction: the subdivision multiplies the cell count" * doctest::test_suite("io")) {
  // what `refinement = 1` does: every element is written as four subcells
  instance::geometry::WriterConfig config;
  config.order = 0;
  config.format = instance::geometry::WriterFormat::Vtk;
  instance::geometry::GeometryWriter geometry(
      "volume", 2, instance::geometry::Shape::Tetrahedron, config, 4);

  std::vector<std::size_t> seenCells;
  std::vector<std::size_t> seenSubcells;
  geometry.addGeometryOutput<double>(
      "v1", {}, false, [&](double* target, std::size_t cell, std::size_t subcell) {
        seenCells.push_back(cell);
        seenSubcells.push_back(subcell);
        target[0] = 0;
      });

  const unit_test::io::TempDir dir;
  auto plan = geometry.makeWriter()(dir.prefix(), 0, 0.0);
  CHECK(inlineValue(plan, "/VTKHDF/NumberOfCells") == 8);

  // the writer function is only invoked once the buffers are materialised
  unit_test::io::runPlan(plan, MPI_COMM_SELF);

  REQUIRE(seenCells.size() == 8);
  CHECK(*std::max_element(seenCells.begin(), seenCells.end()) == 1);
  CHECK(*std::max_element(seenSubcells.begin(), seenSubcells.end()) == 3);
}

} // namespace seissol::unit_test
