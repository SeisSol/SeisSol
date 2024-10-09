/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2017, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Main C++ SeisSol file
 */

#ifndef SEISSOL_H
#define SEISSOL_H

#include <IO/Manager.h>
#include <memory>
#include <string>

#include "utils/logger.h"

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/TimeStepping/LtsLayout.h"
#include "Initializer/Typedefs.h"
#include "Monitoring/FlopCounter.h"
#include "Parallel/Pin.h"
#include "Physics/InstantaneousTimeMirrorManager.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/AsyncIO.h"
#include "ResultWriter/EnergyOutput.h"
#include "ResultWriter/PostProcessor.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/Simulator.h"
#include "Solver/time_stepping/TimeManager.h"
#include "SourceTerm/Manager.h"

namespace seissol {

namespace geometry {
class MeshReader;
} // namespace geometry

/**
 * @todo Initialize rank
 */
class SeisSol {
  public:
  /**
   * Cleanup data structures
   */
  virtual ~SeisSol() { delete m_meshReader; }

  const parallel::Pinning& getPinning() { return pinning; }

  /**
   * Initialize C++ part of the program
   */
  bool init(int argc, char* argv[]);

  /**
   * Finalize SeisSol
   */
  void finalize();

  void loadCheckpoint(const std::string& file);

  initializer::time_stepping::LtsLayout& getLtsLayout() { return m_ltsLayout; }

  initializer::MemoryManager& getMemoryManager() { return *(m_memoryManager.get()); }

  time_stepping::TimeManager& timeManager() { return m_timeManager; }

  Simulator& simulator() { return m_simulator; }

  sourceterm::Manager& sourceTermManager() { return m_sourceTermManager; }

  solver::FreeSurfaceIntegrator& freeSurfaceIntegrator() { return m_freeSurfaceIntegrator; }

  writer::AnalysisWriter& analysisWriter() { return m_analysisWriter; }

  /** Get the post processor module
   */
  writer::PostProcessor& postProcessor() { return m_postProcessor; }

  io::AsyncIO& asyncIO() { return m_asyncIO; }

  /**
   * Get the receiver writer module
   */
  writer::ReceiverWriter& receiverWriter() { return m_receiverWriter; }

  /**
   * Get the energy writer module
   */
  writer::EnergyOutput& energyOutput() { return m_energyOutput; }

  /**
   * Get the flop counter
   */
  monitoring::FlopCounter& flopCounter() { return m_flopCounter; }

  const std::optional<std::string>& getCheckpointLoadFile() { return checkpointLoadFile; }
  /**
   * Reference for timeMirrorManagers to be accessed externally when required
   */
  std::pair<seissol::ITM::InstantaneousTimeMirrorManager,
            seissol::ITM::InstantaneousTimeMirrorManager>&
      getTimeMirrorManagers() {
    return timeMirrorManagers;
  }

  /**
   * Set the mesh reader
   */
  void setMeshReader(seissol::geometry::MeshReader* meshReader) {
    if (m_meshReader != nullptr) {
      logError() << "Mesh reader already initialized";
    }

    m_meshReader = meshReader;
  }

  /**
   * Delete the mesh reader to free memory resources.
   *
   * Should be called after initialization
   */
  void freeMeshReader() {
    delete m_meshReader;
    m_meshReader = 0L;
  }

  /**
   * Get the mesh reader
   */
  const seissol::geometry::MeshReader& meshReader() const { return *m_meshReader; }

  /**
   * Get the mesh reader
   */
  seissol::geometry::MeshReader& meshReader() { return *m_meshReader; }

  seissol::initializer::parameters::SeisSolParameters& getSeisSolParameters() {
    return m_seissolParameters;
  }

  /**
   * Deletes memoryManager. MemoryManager desctructor will destroy LTS Tree and
   * memoryAllocator i.e., the main components of SeisSol. Therefore, call this function
   * at the very end of a program execution
   */
  void deleteMemoryManager() { m_memoryManager.reset(nullptr); }

  GravitationSetup& getGravitationSetup() { return gravitationSetup; }

  /*
   * sets a time stamp for backuping
   * */
  void setBackupTimeStamp(const std::string& stamp);

  /*
   * returns the backup time stamp
   * */
  const std::string& getBackupTimeStamp() { return m_backupTimeStamp; }

  seissol::io::OutputManager& getOutputManager() { return outputManager; }

  private:
  // Note: This HAS to be the first member so that it is initialized before all others!
  // Otherwise it will NOT work.
  // The reason for this is simple yet weird:
  // MPI sets the affinity mask for the process
  // After the first OpenMP call, the OMP runtime sets the pining specified in e.g. OMP_PLACES
  // => Initialize it first, to avoid this.
  parallel::Pinning pinning;

  seissol::io::OutputManager outputManager;

  //! Collection of Parameters
  seissol::initializer::parameters::SeisSolParameters& m_seissolParameters;

  //! Gravitation setup for tsunami boundary condition
  GravitationSetup gravitationSetup;

  //! Async I/O handler (needs to be initialize before other I/O modules)
  io::AsyncIO m_asyncIO;

  //! Mesh Reader
  seissol::geometry::MeshReader* m_meshReader;

  //! Lts Layout
  initializer::time_stepping::LtsLayout m_ltsLayout;

  //! Memory Manager
  std::unique_ptr<initializer::MemoryManager> m_memoryManager{nullptr};

  //! Time Manager
  time_stepping::TimeManager m_timeManager;

  //! Simulator
  Simulator m_simulator;

  //! Source term module
  sourceterm::Manager m_sourceTermManager;

  //! PostProcessor module
  writer::PostProcessor m_postProcessor;

  //! Free surface integrator module
  solver::FreeSurfaceIntegrator m_freeSurfaceIntegrator;

  //! Analysis writer module
  writer::AnalysisWriter m_analysisWriter;

  //! Receiver writer module
  writer::ReceiverWriter m_receiverWriter;

  //! Energy writer module
  writer::EnergyOutput m_energyOutput;

  //! Flop Counter
  monitoring::FlopCounter m_flopCounter;

  //! TimeMirror Managers
  std::pair<seissol::ITM::InstantaneousTimeMirrorManager,
            seissol::ITM::InstantaneousTimeMirrorManager>
      timeMirrorManagers;

  //! time stamp which can be used for backuping files of previous runs
  std::string m_backupTimeStamp{};

  std::optional<std::string> checkpointLoadFile;

  public:
  SeisSol(initializer::parameters::SeisSolParameters& parameters)
      : pinning(), outputManager(*this), m_seissolParameters(parameters), m_meshReader(nullptr),
        m_ltsLayout(parameters),
        m_memoryManager(std::make_unique<initializer::MemoryManager>(*this)), m_timeManager(*this),
        m_analysisWriter(*this), m_receiverWriter(*this), m_energyOutput(*this),
        timeMirrorManagers(*this, *this) {}
};

} // namespace seissol

#endif // SEISSOL_H
