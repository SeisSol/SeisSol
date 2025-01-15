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
#include <Parallel/Pin.h>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace seissol::io::writer::module {

struct BufferPointer {
  size_t size;
  int id;
};

class WriterModule : public seissol::Module, private AsyncWriterModule {
  public:
  WriterModule(const std::string& prefix,
               const ScheduledWriter& settings,
               const parallel::Pinning& pinning);
  void startup();
  void setUp() override;
  void simulationStart() override;
  void syncPoint(double time) override;
  void simulationEnd() override;
  void shutdown() override;

  private:
  int rank;
  std::string prefix;
  unsigned planId;
  AsyncWriter executor;
  std::unordered_map<const void*, BufferPointer> pointerMap;
  std::unordered_map<std::size_t, std::vector<int>> bufferMap;
  ScheduledWriter settings;
  double lastWrite{-1};
  const parallel::Pinning& pinning;
};

} // namespace seissol::io::writer::module

#endif // SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_
