/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#ifndef FAULTWRITER_H
#define FAULTWRITER_H

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

    const int rank = seissol::MPI::mpi.rank();

    wait();

    logInfo(rank) << "Writing faultoutput at time" << utils::nospace << time << ".";

    FaultParam param;
    param.time = time;

    for (unsigned int i = 0; i < m_numVariables; i++) {
      sendBuffer(FaultWriterExecutor::Variables0 + i);
    }

    call(param);

    // Update the timestep count
    m_timestep++;

    m_stopwatch.pause();

    logInfo(rank) << "Writing faultoutput at time" << utils::nospace << time << ". Done.";
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

#endif // FAULTWRITER_H
