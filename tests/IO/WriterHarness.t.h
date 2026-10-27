// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_IO_WRITERHARNESS_T_H_
#define SEISSOL_TESTS_IO_WRITERHARNESS_T_H_

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Instruction.h"
#include "IO/Writer/Module/AsyncWriter.h"
#include "IO/Writer/Writer.h"

#include <algorithm>
#include <async/ExecInfo.h>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <mpi.h>
#include <set>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace seissol::unit_test::io {

/**
 * A stand-in for the ASYNC executor, so that a write plan can be run inside a unit test.
 *
 * It performs the same three steps WriterModule does: materialize every managed buffer, hand out
 * ids for them, and expose them again through ExecInfo. The plan is serialized and re-parsed in
 * between, exactly as it is when it travels to the executor, so the round trip through YAML is
 * part of what gets exercised here.
 */
class LocalExecInfo : public async::ExecInfo {
  public:
  //! @param comm The ranks that execute, as ASYNC hands them to the executor
  explicit LocalExecInfo(MPI_Comm comm = MPI_COMM_SELF) { setCommInternal(comm); }

  [[nodiscard]] const void* buffer(unsigned int id) const override { return buffers_[id].data(); }

  int addBuffer(std::size_t size) {
    buffers_.emplace_back(size);
    addBufferInternal(size);
    return static_cast<int>(buffers_.size()) - 1;
  }

  void resizeBuffer(int id, std::size_t size) {
    buffers_[id].resize(size);
    resizeBufferInternal(id, size);
  }

  std::vector<char>& bufferData(int id) { return buffers_[id]; }

  /**
   * Let a buffer stand for the contributions of several ranks, as the executor of ASYNC mode NODE
   * holds them: one after another, in @p parts bytes each
   */
  void setParts(int id, std::vector<std::size_t> parts) { parts_[id] = std::move(parts); }

  [[nodiscard]] std::vector<std::size_t> bufferParts(unsigned int id) const override {
    const auto it = parts_.find(static_cast<int>(id));
    return it == parts_.end() ? async::ExecInfo::bufferParts(id) : it->second;
  }

  private:
  std::vector<std::vector<char>> buffers_;
  std::map<int, std::vector<std::size_t>> parts_;
};

/**
 * Materializes the managed buffers of @p writer , serializes the plan, parses it back and executes
 * it on @p comm . Returns the serialized plan, so that a test can inspect it as well.
 *
 * With @p runFiles , the parsed plan is carried out directly against that record of the files a
 * run has written, so that a test can play several runs one after another; otherwise it goes
 * through the executor, which keeps that record for the whole process.
 */
inline std::string runPlan(seissol::io::writer::Writer& writer,
                           MPI_Comm comm,
                           seissol::io::writer::file::RunFiles* runFiles = nullptr) {
  using namespace seissol::io::writer;

  LocalExecInfo info(comm);
  // buffer 0 is where WriterModule puts the serialized plan; the executor reads it from there
  const auto planId = info.addBuffer(0);

  std::set<DataSource*> handled;
  for (const auto& instruction : writer.getInstructions()) {
    for (const auto& source : instruction->dataSources()) {
      if (!source->distributed() || handled.count(source.get()) > 0) {
        continue;
      }
      handled.emplace(source.get());
      if (auto* adhoc = dynamic_cast<AdhocBuffer*>(source.get()); adhoc != nullptr) {
        const auto id = info.addBuffer(adhoc->getTargetSize());
        adhoc->setData(info.bufferData(id).data());
        source->assignId(id);
      } else {
        const auto id = info.addBuffer(source->getLocalSize());
        auto& target = info.bufferData(id);
        const auto* local = static_cast<const char*>(source->getLocalPointer());
        std::copy(local, local + target.size(), target.begin());
        source->assignId(id);
      }
    }
  }

  // ids are assigned by now, so the plan is final
  const auto plan = writer.serialize();
  info.resizeBuffer(planId, plan.size());
  std::ranges::copy(plan, info.bufferData(planId).begin());

  if (runFiles != nullptr) {
    WriteInstance instance(comm, runFiles);
    Writer(plan).write(info, instance);
    instance.close();
    return plan;
  }

  // hand it to the real executor rather than writing the plan directly, so that the plan travels
  // the way it does in production
  module::AsyncWriter executor;
  executor.execInit(info, module::AsyncWriterInit{});
  executor.exec(info, module::AsyncWriterExec{});
  executor.execWait(info);
  executor.closeFiles();

  return plan;
}

//! A temporary directory that removes itself again.
struct TempDir {
  std::string path;

  TempDir() {
    // mkdtemp rewrites its argument, so it needs a fresh buffer on every call
    std::string buffer = "/tmp/seissoliotestXXXXXX";
    const char* created = mkdtemp(buffer.data());
    if (created == nullptr) {
      throw std::runtime_error("could not create a temporary directory");
    }
    path = created;
  }
  ~TempDir() {
    std::error_code error;
    std::filesystem::remove_all(path, error);
  }
  TempDir(const TempDir&) = delete;
  TempDir(TempDir&&) = delete;
  auto operator=(const TempDir&) -> TempDir& = delete;
  auto operator=(TempDir&&) -> TempDir& = delete;

  [[nodiscard]] std::string prefix() const { return path + "/out"; }
};

/**
 * Reads a table of samples, each a compound of the doubles @p quantities , as the flat array of
 * numbers it is in memory: (sample, point, quantity). HDF5 converts a compound only into another
 * compound, matching the members by name, so the target type is spelled out here.
 */
inline std::vector<double> readSampleTable(seissol::io::reader::file::Hdf5Reader& hdf5,
                                           const std::string& name,
                                           const std::vector<std::string>& quantities) {
  using namespace seissol::io::datatype;
  std::vector<StructDatatype::MemberInfo> members(quantities.size());
  for (std::size_t i = 0; i < quantities.size(); ++i) {
    members[i] = StructDatatype::MemberInfo{
        .name = quantities[i], .offset = i * sizeof(double), .datatype = inferDatatype<double>()};
  }
  const auto rows = hdf5.dataCount(name);
  std::vector<double> values(rows * hdf5.dataRowSize(name) * quantities.size());
  hdf5.readDataRaw(values.data(), name, rows, std::make_shared<StructDatatype>(members));
  return values;
}

/**
 * Two tetrahedra sharing a face, in the vertex order the output writes them: the point projector
 * is handed one cell at a time and fills the four corners.
 */
constexpr double TestVertices[2][4][3] = {
    {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
    {{1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.0}}};

} // namespace seissol::unit_test::io

#endif // SEISSOL_TESTS_IO_WRITERHARNESS_T_H_
