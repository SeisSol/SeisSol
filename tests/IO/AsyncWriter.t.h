// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Module/AsyncWriter.h"
#include "IO/Writer/Writer.h"
#include "WriterHarness.t.h"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using namespace seissol::io;
using namespace seissol::io::writer;

//! A plan writing one managed buffer, plus the values that buffer is expected to hold.
std::pair<Writer, std::vector<std::int64_t>>
    makePlan(const std::string& file, const std::string& name, std::int64_t base) {
  std::vector<std::int64_t> expected(16);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    expected[i] = base + static_cast<std::int64_t>(i);
  }

  auto source = GeneratedBuffer::createElementwise<std::int64_t>(
      expected.size(), 1, {}, [base](std::int64_t* target, std::size_t index) {
        target[0] = base + static_cast<std::int64_t>(index);
      });

  Writer writer;
  writer.addInstruction(
      std::make_shared<instructions::Hdf5DataWrite>(instructions::Hdf5Location(file, {"group"}),
                                                    name,
                                                    source,
                                                    datatype::inferDatatype<std::int64_t>()));
  return {std::move(writer), expected};
}

} // namespace

// ---------------------------------------------------------------------------
// The executor side: everything reaches it through ExecInfo buffers
// ---------------------------------------------------------------------------

/**
 * WriterModule and AsyncWriter only share ExecInfo: buffer 0 carries the serialised plan, the
 * remaining buffers carry the data, and the plan refers to them by the ids WriterModule assigned.
 * If those two ends ever disagree about an id, the executor writes whatever happens to sit at that
 * index -- silently, since every buffer is just bytes.
 */
TEST_CASE("IO/AsyncWriter: the plan and its buffers travel through ExecInfo" *
          doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;
  const auto file = dir.path + "/async.h5";

  auto [writer, expected] = makePlan(file, "values", 100);

  io::LocalExecInfo info;
  const auto planId = info.addBuffer(0);
  CHECK(planId == 0);

  // materialise the managed buffer and tell the plan where it went, as WriterModule does
  auto sources = writer.getInstructions().front()->dataSources();
  REQUIRE(sources.size() == 1);
  auto* adhoc = dynamic_cast<AdhocBuffer*>(sources.front().get());
  REQUIRE(adhoc != nullptr);
  const auto dataId = info.addBuffer(adhoc->getTargetSize());
  CHECK(dataId == 1);
  adhoc->setData(info.bufferData(dataId).data());
  sources.front()->assignId(dataId);

  const auto plan = writer.serialize();
  info.resizeBuffer(planId, plan.size());
  std::copy(plan.begin(), plan.end(), info.bufferData(planId).begin());

  // the id has to survive into the serialised plan, otherwise the executor cannot find the data
  CHECK(plan.find("id: " + std::to_string(dataId)) != std::string::npos);

  module::AsyncWriter executor;
  executor.setComm(MPI_COMM_SELF);
  executor.execInit(info, module::AsyncWriterInit{});
  executor.exec(info, module::AsyncWriterExec{});
  executor.execWait(info);
  executor.finalize();

  REQUIRE(std::filesystem::exists(file));
  reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
  hdf5.openFile(file);
  hdf5.openGroup("group");
  const auto values = hdf5.readData<std::int64_t>("values");
  REQUIRE(values.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(values[i] == expected[i]);
  }
  hdf5.closeGroup();
  hdf5.closeFile();
}

/**
 * The executor keeps its Writer and its WriteInstance as members and is reused across output
 * steps, so a leftover from one plan must not leak into the next one.
 */
TEST_CASE("IO/AsyncWriter: one executor runs several plans" * doctest::test_suite("io")) {
  const unit_test::io::TempDir dir;

  module::AsyncWriter executor;
  executor.setComm(MPI_COMM_SELF);

  for (std::int64_t step = 0; step < 3; ++step) {
    const auto file = dir.path + "/async-" + std::to_string(step) + ".h5";
    auto [writer, expected] = makePlan(file, "values", 1000 * (step + 1));

    io::LocalExecInfo info;
    const auto planId = info.addBuffer(0);
    auto sources = writer.getInstructions().front()->dataSources();
    auto* adhoc = dynamic_cast<AdhocBuffer*>(sources.front().get());
    const auto dataId = info.addBuffer(adhoc->getTargetSize());
    adhoc->setData(info.bufferData(dataId).data());
    sources.front()->assignId(dataId);

    const auto plan = writer.serialize();
    info.resizeBuffer(planId, plan.size());
    std::copy(plan.begin(), plan.end(), info.bufferData(planId).begin());

    executor.execInit(info, module::AsyncWriterInit{});
    executor.exec(info, module::AsyncWriterExec{});
    executor.execWait(info);

    reader::file::Hdf5Reader hdf5(MPI_COMM_SELF);
    hdf5.openFile(file);
    hdf5.openGroup("group");
    const auto values = hdf5.readData<std::int64_t>("values");
    REQUIRE(values.size() == expected.size());
    CHECK(values.front() == expected.front());
    CHECK(values.back() == expected.back());
    hdf5.closeGroup();
    hdf5.closeFile();
  }

  executor.finalize();
}

/**
 * The plan is written into a buffer of exactly its own length, so a serialisation that is not
 * stable under a round trip would truncate or corrupt it on the next step.
 */
TEST_CASE("IO/AsyncWriter: the serialised plan is stable" * doctest::test_suite("io")) {
  auto [writer, expected] = makePlan("some.h5", "values", 0);
  writer.getInstructions().front()->dataSources().front()->assignId(7);

  const auto once = writer.serialize();
  Writer restored(once);
  const auto twice = restored.serialize();
  CHECK(once == twice);

  Writer restoredAgain(twice);
  CHECK(restoredAgain.getInstructions().size() == writer.getInstructions().size());
  CHECK(restoredAgain.serialize() == once);
}

} // namespace seissol::unit_test
