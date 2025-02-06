// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_FAULTWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_FAULTWRITER_H_

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include "utils/logger.h"

#include "async/Module.h"

#include "FaultWriterExecutor.h"
#include "Modules/Module.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/Stopwatch.h"

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
  bool m_enabled{false};

  /** The asynchronous executor */
  FaultWriterExecutor m_executor;

  /** Total number of variables */
  unsigned int m_numVariables{0};

  /** The current output time step */
  unsigned int m_timestep{0};

  /** Frontend stopwatch */
  Stopwatch m_stopwatch;

  dr::output::OutputManager* callbackObject{nullptr};

  public:
  FaultWriter(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance)

  {}

  /**
   * Called by ASYNC on all ranks
   */
  void setUp() override;

  void setTimestep(unsigned int timestep) { m_timestep = timestep; }

  void init(const unsigned int* cells,
            const double* vertices,
            const unsigned int* faultTags,
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
  [[nodiscard]] unsigned int timestep() const { return m_timestep; }

  void write(double time) {
    SCOREP_USER_REGION("FaultWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION)

    if (!m_enabled) {
      logError() << "Trying to write fault output, but fault output is not enabled";
    }

    m_stopwatch.start();

    wait();

    logInfo() << "Writing faultoutput at time" << utils::nospace << time << ".";

    FaultParam param;
    param.time = time;

    for (unsigned int i = 0; i < m_numVariables; i++) {
      sendBuffer(FaultWriterExecutor::Variables0 + i);
    }

    call(param);

    // Update the timestep count
    m_timestep++;

    m_stopwatch.pause();

    logInfo() << "Writing faultoutput at time" << utils::nospace << time << ". Done.";
  }

  void close() {
    if (m_enabled) {
      wait();
    }

    finalize();

    if (!m_enabled) {
      return;
    }

    m_stopwatch.printTime("Time fault writer frontend:");
  }

  void tearDown() override { m_executor.finalize(); }

  void setupCallbackObject(dr::output::OutputManager* faultOutputManager) {
    callbackObject = faultOutputManager;
  }

  //
  // Hooks
  //
  void simulationStart() override;

  void syncPoint(double currentTime) override;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_FAULTWRITER_H_
