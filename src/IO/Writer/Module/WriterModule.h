// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_
#define SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_

#include "AsyncWriter.h"
#include "Modules/Module.h"
#include "Parallel/Pin.h"

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::io::writer::module {

struct BufferPointer {
  size_t size{0};
  int id{-1};
};

class WriterModule : public seissol::Module, private AsyncWriterModule {
  public:
  WriterModule(const std::string& prefix,
               const ScheduledWriter& settings,
               const parallel::Pinning& pinning,
               SeisSol& seissolInstance);
  void startup();
  void simulationStart(std::optional<double> checkpointTime) override;
  void syncPoint(double time) override;
  void simulationEnd() override;
  void shutdown() override;

  private:
  void setUp() override;

  int rank;
  std::string prefix;
  unsigned planId{std::numeric_limits<unsigned>::max()};
  AsyncWriter executor;
  std::unordered_map<const void*, BufferPointer> pointerMap;
  std::unordered_map<std::size_t, std::vector<int>> bufferMap;
  ScheduledWriter settings;
  double lastWrite{-1};
  const parallel::Pinning& pinning;

  // TODO: remove?
  SeisSol& seissolInstance;
};

} // namespace seissol::io::writer::module

#endif // SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_
