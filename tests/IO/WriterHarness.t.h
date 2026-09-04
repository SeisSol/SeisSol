// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_IO_WRITERHARNESS_T_H_
#define SEISSOL_TESTS_IO_WRITERHARNESS_T_H_

#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Instruction.h"
#include "IO/Writer/Module/AsyncWriter.h"
#include "IO/Writer/Writer.h"

#include <async/ExecInfo.h>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <mpi.h>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test::io {

/**
 * A stand-in for the ASYNC executor, so that a write plan can be run inside a unit test.
 *
 * It performs the same three steps WriterModule does: materialise every managed buffer, hand out
 * ids for them, and expose them again through ExecInfo. The plan is serialised and re-parsed in
 * between, exactly as it is when it travels to the executor, so the round trip through YAML is
 * part of what gets exercised here.
 */
class LocalExecInfo : public async::ExecInfo {
  public:
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

  private:
  std::vector<std::vector<char>> buffers_;
};

/**
 * Materialises the managed buffers of @p writer , serialises the plan, parses it back and executes
 * it on @p comm . Returns the serialised plan, so that a test can inspect it as well.
 */
inline std::string runPlan(seissol::io::writer::Writer& writer, MPI_Comm comm) {
  using namespace seissol::io::writer;

  LocalExecInfo info;
  // buffer 0 is where WriterModule puts the serialised plan; the executor reads it from there
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
  std::copy(plan.begin(), plan.end(), info.bufferData(planId).begin());

  // hand it to the real executor rather than calling beginWrite directly, so that the plan
  // travels the way it does in production
  module::AsyncWriter executor;
  executor.setComm(comm);
  executor.execInit(info, module::AsyncWriterInit{});
  executor.exec(info, module::AsyncWriterExec{});
  executor.execWait(info);
  executor.finalize();

  writer.endWrite();
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
 * Two tetrahedra sharing a face, in the vertex order the output writes them: the point projector
 * is handed one cell at a time and fills the four corners.
 */
constexpr double TestVertices[2][4][3] = {
    {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
    {{1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 0.0}}};

} // namespace seissol::unit_test::io

#endif // SEISSOL_TESTS_IO_WRITERHARNESS_T_H_
