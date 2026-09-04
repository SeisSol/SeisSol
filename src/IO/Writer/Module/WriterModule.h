// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_
#define SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_

#include "AsyncWriter.h"
#include "BufferRegistry.h"
#include "Modules/Module.h"
#include "Parallel/Pin.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <string>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::io::writer::module {

/**
 * @brief The number of outputs a run at @p startTime has behind it, i.e. the counter its first
 * output should carry.
 *
 * Outputs are numbered from zero and happen every @p interval, so a run resuming from a
 * checkpoint has to pick up where the previous one stopped instead of starting over and
 * overwriting its files. The comparison is made with a tolerance, since a checkpoint time that
 * should sit exactly on an output can arrive a rounding error below it.
 */
inline std::size_t outputCountBefore(double startTime, double interval) {
  if (!(startTime > 0) || !(interval > 0)) {
    return 0;
  }
  constexpr double Tolerance = 1e-9;
  const auto elapsed = startTime / interval;
  // the output at startTime itself has already happened, hence the + 1
  return static_cast<std::size_t>(std::floor(elapsed + Tolerance * std::max(1.0, elapsed))) + 1;
}

class WriterModule : public seissol::Module, private AsyncWriterModule, private BufferAllocator {
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

  // BufferAllocator, forwarding to ASYNC
  int addBuffer(const void* pointer, std::size_t size) override;
  void resizeBuffer(int id, const void* pointer, std::size_t size) override;
  void* managedBuffer(int id) override;

  int rank_;
  std::string prefix_;
  unsigned planId_{std::numeric_limits<unsigned>::max()};
  AsyncWriter executor_;
  BufferRegistry registry_{*this};
  ScheduledWriter settings_;
  double lastWrite_{-1};
  std::size_t writeCount_{0};
  const parallel::Pinning& pinning_;

  // TODO: remove?
  SeisSol& seissolInstance_;
};

} // namespace seissol::io::writer::module

#endif // SEISSOL_SRC_IO_WRITER_MODULE_WRITERMODULE_H_
