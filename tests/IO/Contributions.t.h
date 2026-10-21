// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/Writer/File/RunFiles.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Dimension.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "WriterHarness.t.h"

#include <cstddef>
#include <cstring>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace contributionstest {
using namespace seissol::io;
using namespace seissol::io::writer;

//! The value at a point (along the distributed dimension) and an entry of the other one.
inline double entryValue(std::size_t point, std::size_t entry) {
  return (100.0 * static_cast<double>(point)) + static_cast<double>(entry);
}

/**
 * Writes @p data, of the shape @p dimensions, as the executor of ASYNC mode NODE does: from one
 * buffer that holds the contributions of several ranks, of @p parts bytes each. Returns what the
 * file holds afterwards.
 */
inline std::vector<double> writeContributions(const std::vector<double>& data,
                                              const std::vector<std::size_t>& parts,
                                              const std::vector<Dimension>& dimensions) {
  const unit_test::io::TempDir dir;
  const auto file = dir.prefix() + ".h5";

  // the plan, as it reaches the executor: the data lives in buffer 0
  auto source = std::make_shared<WriteBuffer>(
      data.data(), data.size(), datatype::inferDatatype<double>(), dimensions);
  source->assignId(0);
  Writer writer;
  writer.addInstruction(
      std::make_shared<instructions::Hdf5DataWrite>(instructions::Hdf5Location(file, {"table"}),
                                                    "values",
                                                    source,
                                                    datatype::inferDatatype<double>()));
  const auto plan = writer.serialize();

  unit_test::io::LocalExecInfo info(MPI_COMM_SELF);
  const auto id = info.addBuffer(data.size() * sizeof(double));
  std::memcpy(info.bufferData(id).data(), data.data(), data.size() * sizeof(double));
  info.setParts(id, parts);

  file::RunFiles runFiles;
  WriteInstance instance(MPI_COMM_SELF, &runFiles);
  Writer(plan).write(info, instance);
  instance.close();

  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(file);
  hdf5.openGroup("table");
  auto values = hdf5.readData<double>("values");
  hdf5.closeGroup();
  hdf5.closeFile();
  return values;
}

} // namespace contributionstest

using namespace contributionstest;

TEST_CASE("IO/Hdf5: the contributions of the ranks of a group are put together along a "
          "distributed dimension that is not the leading one" *
          doctest::test_suite("io")) {
  // a table of entries (samples, say) by points, from four ranks, one of them without points: each
  // rank contributes a block of all entries by its points
  constexpr std::size_t Entries = 3;
  const std::vector<std::size_t> points{2, 1, 0, 3};

  std::vector<double> data;
  std::vector<std::size_t> parts;
  std::size_t total = 0;
  for (const auto rankPoints : points) {
    for (std::size_t entry = 0; entry < Entries; ++entry) {
      for (std::size_t point = 0; point < rankPoints; ++point) {
        data.push_back(entryValue(total + point, entry));
      }
    }
    parts.push_back(Entries * rankPoints * sizeof(double));
    total += rankPoints;
  }

  const auto values =
      writeContributions(data, parts, {Dimension::replicated(Entries), Dimension::distributed()});

  REQUIRE(values.size() == Entries * total);
  for (std::size_t entry = 0; entry < Entries; ++entry) {
    for (std::size_t point = 0; point < total; ++point) {
      CHECK(values[(entry * total) + point] == doctest::Approx(entryValue(point, entry)));
    }
  }
}

TEST_CASE("IO/Hdf5: the contributions along a leading distributed dimension are the block" *
          doctest::test_suite("io")) {
  // points by entries: the blocks of the ranks follow each other as rows
  constexpr std::size_t Entries = 2;
  const std::vector<std::size_t> points{1, 3};

  std::vector<double> data;
  std::vector<std::size_t> parts;
  std::size_t total = 0;
  for (const auto rankPoints : points) {
    for (std::size_t point = 0; point < rankPoints; ++point) {
      for (std::size_t entry = 0; entry < Entries; ++entry) {
        data.push_back(entryValue(total + point, entry));
      }
    }
    parts.push_back(rankPoints * Entries * sizeof(double));
    total += rankPoints;
  }

  const auto values =
      writeContributions(data, parts, {Dimension::distributed(), Dimension::replicated(Entries)});

  REQUIRE(values.size() == total * Entries);
  for (std::size_t point = 0; point < total; ++point) {
    for (std::size_t entry = 0; entry < Entries; ++entry) {
      CHECK(values[(point * Entries) + entry] == doctest::Approx(entryValue(point, entry)));
    }
  }
}

} // namespace seissol::unit_test
