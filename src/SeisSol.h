// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_SEISSOL_H_
#define SEISSOL_SRC_SEISSOL_H_

#include <Common/Executor.h>
#include <IO/Manager.h>
#include <memory>
#include <optional>
#include <string>

#include "utils/env.h"
#include "utils/logger.h"

#include "Initializer/Parameters/SeisSolParameters.h"
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
#include "Solver/TimeStepping/TimeManager.h"
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

  Executor executionPlace(std::size_t clusterSize) {
    constexpr auto DefaultDevice = isDeviceOn() ? Executor::Device : Executor::Host;
    if (executionPlaceCutoff.has_value()) {
      if (executionPlaceCutoff.value() <= clusterSize) {
        return DefaultDevice;
      } else {
        return Executor::Host;
      }
    } else {
      return DefaultDevice;
    }
  }

  void setExecutionPlaceCutoff(std::size_t size);

  initializer::MemoryManager& getMemoryManager() { return *m_memoryManager; }

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
    m_meshReader = nullptr;
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

  const seissol::initializer::parameters::SeisSolParameters& getSeisSolParameters() const {
    return m_seissolParameters;
  }

  /**
   * Deletes memoryManager. MemoryManager desctructor will destroy all storage structures and
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

  utils::Env& env() { return m_env; }

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
  seissol::geometry::MeshReader* m_meshReader{nullptr};

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
  std::string m_backupTimeStamp;

  std::optional<std::string> checkpointLoadFile;

  std::optional<std::size_t> executionPlaceCutoff;

  utils::Env m_env;

  public:
  SeisSol(initializer::parameters::SeisSolParameters& parameters, const utils::Env& env)
      : outputManager(*this), m_seissolParameters(parameters),
        m_memoryManager(std::make_unique<initializer::MemoryManager>(*this)), m_timeManager(*this),
        m_analysisWriter(*this), m_receiverWriter(*this), m_energyOutput(*this),
        timeMirrorManagers(*this, *this), m_env(env) {}

  SeisSol(const SeisSol&) = delete;
  SeisSol(SeisSol&&) = delete;
  auto operator=(const SeisSol&) = delete;
  auto operator=(SeisSol&&) = delete;
};

} // namespace seissol

#endif // SEISSOL_SRC_SEISSOL_H_
