// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_FAULTWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_FAULTWRITER_H_

#include "FaultWriterExecutor.h"
#include "Modules/Module.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <async/Module.h>
#include <utils/logger.h>

namespace seissol {
class SeisSol;
namespace dr::output {
class OutputManager;
} // namespace dr::output
} // namespace seissol

namespace seissol::writer {

class FaultWriter : private async::Module<FaultWriterExecutor, FaultInitParam, FaultParam>,
                    public seissol::Module {
  private:
  seissol::SeisSol& seissolInstance;

  /** Is enabled? */
  bool enabled_{false};

  /** The asynchronous executor */
  FaultWriterExecutor executor_;

  /** Total number of variables */
  unsigned int numVariables_{0};

  /** The current output time step */
  unsigned int timestep_{0};

  /** Frontend stopwatch */
  Stopwatch stopwatch_;

  dr::output::OutputManager* callbackObject{nullptr};

  public:
  explicit FaultWriter(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance)

  {}

  /**
   * Called by ASYNC on all ranks
   */
  void setUp() override;

  void setTimestep(unsigned int timestep) { timestep_ = timestep; }

  void init(const unsigned int* cells,
            const double* vertices,
            const unsigned int* faultTags,
            const unsigned int* ids,
            unsigned int nCells,
            unsigned int nVertices,
            const int* outputMask,
            const real** dataBuffer,
            const char* outputPrefix,
            double interval,
            xdmfwriter::BackendType backend,
            const std::string& backupTimeStamp);

  /**
   * @return The current time step of the fault output
   */
  [[nodiscard]] unsigned int timestep() const { return timestep_; }

  void write(double time) {
    SCOREP_USER_REGION("FaultWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION)

    if (!enabled_) {
      logError() << "Trying to write fault output, but fault output is not enabled";
    }

    stopwatch_.start();

    wait();

    logInfo() << "Writing faultoutput at time" << utils::nospace << time << ".";

    FaultParam param;
    param.time = time;

    for (unsigned int i = 0; i < numVariables_; i++) {
      sendBuffer(FaultWriterExecutor::Variables0 + i);
    }

    call(param);

    // Update the timestep count
    timestep_++;

    stopwatch_.pause();

    logInfo() << "Writing faultoutput at time" << utils::nospace << time << ". Done.";
  }

  void close() {
    if (enabled_) {
      wait();
    }

    finalize();

    if (!enabled_) {
      return;
    }

    stopwatch_.printTime("Time fault writer frontend:");
  }

  void tearDown() override { executor_.finalize(); }

  void setupCallbackObject(dr::output::OutputManager* faultOutputManager) {
    callbackObject = faultOutputManager;
  }

  //
  // Hooks
  //
  void simulationStart(std::optional<double> checkpointTime) override;

  void syncPoint(double currentTime) override;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_FAULTWRITER_H_
