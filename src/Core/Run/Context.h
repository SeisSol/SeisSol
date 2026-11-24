// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_RUN_CONTEXT_H_
#define SEISSOL_SRC_CORE_RUN_CONTEXT_H_

#include "Core/IO/AsyncIO.h"
#include "Core/IO/Manager.h"
#include "Core/IO/OutputContext.h"
#include "Core/Modules/Modules.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Parallel/Pin.h"
#include "Core/Parallel/Runtime/HostActor.h"

#include <string>
#include <utility>
#include <utils/env.h>

namespace seissol {

/**
 * The services of a run, which its parts use whatever parts it has: its environment, the pinning of
 * its threads, its hooks, its asynchronous IO and its outputs, its FLOP counter, the actor of its
 * work on the host for the device, and its ID and the prefix of its files. The composition
 * (SeisSol) is one, and holds the parts besides.
 */
class RunContext {
  public:
  explicit RunContext(utils::Env env)
      : env_(std::move(env)), asyncIO_(env_),
        outputManager_(io::OutputContext{.pinning = pinning_,
                                         .modules = modules_,
                                         .asyncIO = asyncIO_,
                                         .env = env_,
                                         .runId = runId_,
                                         .outputPrefix = outputPrefix_}) {}

  // the outputs refer to the services
  RunContext(const RunContext&) = delete;
  RunContext(RunContext&&) = delete;
  auto operator=(const RunContext&) = delete;
  auto operator=(RunContext&&) = delete;
  virtual ~RunContext() = default;

  const parallel::Pinning& pinning() { return pinning_; }

  utils::Env& env() { return env_; }

  Modules& modules() { return modules_; }

  io::AsyncIO& asyncIO() { return asyncIO_; }

  io::OutputManager& outputManager() { return outputManager_; }

  monitoring::FlopCounter& flopCounter() { return flopCounter_; }

  /**
   * Runs the work on the host that reads what the device has written (e.g. the samples of the
   * receivers), for all parts of the run.
   */
  parallel::runtime::HostActor& hostActor() { return hostActor_; }

  /**
   * Sets the ID of this run (see Core/Common/RunId.h); the same on every rank.
   */
  void setRunId(const std::string& runId) { runId_ = runId; }

  /**
   * The ID of this run, which every file of its output carries.
   */
  [[nodiscard]] const std::string& runId() const { return runId_; }

  /**
   * Sets what the names of the output files start with, their directory included (see
   * initprocedure::createOutputDirectory).
   */
  void setOutputPrefix(const std::string& outputPrefix) { outputPrefix_ = outputPrefix; }

  /**
   * What the names of the output files start with, their directory included.
   */
  [[nodiscard]] const std::string& outputPrefix() const { return outputPrefix_; }

  protected:
  // first, so that it outlives the queues of all parts of the run
  parallel::runtime::HostActor hostActor_;
  // the composition sets them up (see SeisSol::init)
  parallel::Pinning pinning_;
  utils::Env env_;
  Modules modules_;
  std::string runId_;
  std::string outputPrefix_;
  io::AsyncIO asyncIO_;
  io::OutputManager outputManager_;
  monitoring::FlopCounter flopCounter_;
};

} // namespace seissol

#endif // SEISSOL_SRC_CORE_RUN_CONTEXT_H_
