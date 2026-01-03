// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_SEISSOL_H_
#define SEISSOL_SRC_SEISSOL_H_

#include "Common/Executor.h"
#include "IO/Manager.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Monitoring/FlopCounter.h"
#include "Parallel/Pin.h"
#include "Physics/InstantaneousTimeMirrorManager.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/AsyncIO.h"
#include "ResultWriter/EnergyOutput.h"
#include "ResultWriter/FaultWriter.h"
#include "ResultWriter/FreeSurfaceWriter.h"
#include "ResultWriter/PostProcessor.h"
#include "ResultWriter/WaveFieldWriter.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/Simulator.h"
#include "Solver/TimeStepping/TimeManager.h"
#include "SourceTerm/Manager.h"

#include <memory>
#include <optional>
#include <string>
#include <utils/env.h>
#include <utils/logger.h>

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
  virtual ~SeisSol() { delete meshReader_; }

  const parallel::Pinning& getPinning() { return pinning_; }

  /**
   * Initialize C++ part of the program
   */
  bool init();

  /**
   * Finalize SeisSol
   */
  void finalize();

  void loadCheckpoint(const std::string& file);

  Executor executionPlace(std::size_t clusterSize) {
    constexpr auto DefaultDevice = isDeviceOn() ? Executor::Device : Executor::Host;
    if (executionPlaceCutoff_.has_value()) {
      if (executionPlaceCutoff_.value() <= clusterSize) {
        return DefaultDevice;
      } else {
        return Executor::Host;
      }
    } else {
      return DefaultDevice;
    }
  }

  void setExecutionPlaceCutoff(std::size_t size);

  initializer::MemoryManager& getMemoryManager() { return *memoryManager_; }

  time_stepping::TimeManager& timeManager() { return timeManager_; }

  Simulator& simulator() { return simulator_; }

  sourceterm::Manager& sourceTermManager() { return sourceTermManager_; }

  solver::FreeSurfaceIntegrator& freeSurfaceIntegrator() { return freeSurfaceIntegrator_; }

  writer::FreeSurfaceWriter& freeSurfaceWriter() { return freeSurfaceWriter_; }

  writer::AnalysisWriter& analysisWriter() { return analysisWriter_; }

  /** Get the post processor module
   */
  writer::PostProcessor& postProcessor() { return postProcessor_; }

  io::AsyncIO& asyncIO() { return asyncIO_; }

  /**
   * Get the wave field writer module
   */
  writer::WaveFieldWriter& waveFieldWriter() { return waveFieldWriter_; }

  /**
   * Get the fault writer module
   */
  writer::FaultWriter& faultWriter() { return faultWriter_; }

  /**
   * Get the receiver writer module
   */
  writer::ReceiverWriter& receiverWriter() { return receiverWriter_; }

  /**
   * Get the energy writer module
   */
  writer::EnergyOutput& energyOutput() { return energyOutput_; }

  /**
   * Get the flop counter
   */
  monitoring::FlopCounter& flopCounter() { return flopCounter_; }

  const std::optional<std::string>& getCheckpointLoadFile() { return checkpointLoadFile_; }
  /**
   * Reference for timeMirrorManagers to be accessed externally when required
   */
  std::pair<seissol::ITM::InstantaneousTimeMirrorManager,
            seissol::ITM::InstantaneousTimeMirrorManager>&
      getTimeMirrorManagers() {
    return timeMirrorManagers_;
  }

  /**
   * Set the mesh reader
   */
  void setMeshReader(seissol::geometry::MeshReader* meshReader) {
    if (meshReader_ != nullptr) {
      logError() << "Mesh reader already initialized";
    }

    meshReader_ = meshReader;
  }

  /**
   * Delete the mesh reader to free memory resources.
   *
   * Should be called after initialization
   */
  void freeMeshReader() {
    delete meshReader_;
    meshReader_ = nullptr;
  }

  /**
   * Get the mesh reader
   */
  const seissol::geometry::MeshReader& meshReader() const { return *meshReader_; }

  /**
   * Get the mesh reader
   */
  seissol::geometry::MeshReader& meshReader() { return *meshReader_; }

  const seissol::initializer::parameters::SeisSolParameters& getSeisSolParameters() const {
    return seissolParameters_;
  }

  /**
   * Deletes memoryManager. MemoryManager desctructor will destroy all storage structures and
   * memoryAllocator i.e., the main components of SeisSol. Therefore, call this function
   * at the very end of a program execution
   */
  void deleteMemoryManager() { memoryManager_.reset(nullptr); }

  GravitationSetup& getGravitationSetup() { return gravitationSetup_; }

  /*
   * sets a time stamp for backuping
   * */
  void setBackupTimeStamp(const std::string& stamp);

  void setTimestepScale(double scale) { timestepScale_ = scale; }

  double getTimestepScale() const { return timestepScale_; }

  /*
   * returns the backup time stamp
   * */
  const std::string& getBackupTimeStamp() { return backupTimeStamp_; }

  seissol::io::OutputManager& getOutputManager() { return outputManager_; }

  utils::Env& env() { return env_; }

  private:
  // Note: This HAS to be the first member so that it is initialized before all others!
  // Otherwise it will NOT work.
  // The reason for this is simple yet weird:
  // MPI sets the affinity mask for the process
  // After the first OpenMP call, the OMP runtime sets the pining specified in e.g. OMP_PLACES
  // => Initialize it first, to avoid this.
  parallel::Pinning pinning_;

  seissol::io::OutputManager outputManager_;

  //! Collection of Parameters
  const seissol::initializer::parameters::SeisSolParameters& seissolParameters_;

  //! Gravitation setup for tsunami boundary condition
  GravitationSetup gravitationSetup_;

  //! Async I/O handler (needs to be initialize before other I/O modules)
  io::AsyncIO asyncIO_;

  //! Mesh Reader
  seissol::geometry::MeshReader* meshReader_{nullptr};

  //! Memory Manager
  std::unique_ptr<initializer::MemoryManager> memoryManager_{nullptr};

  //! Time Manager
  time_stepping::TimeManager timeManager_;

  //! Simulator
  Simulator simulator_;

  //! Source term module
  sourceterm::Manager sourceTermManager_;

  //! PostProcessor module
  writer::PostProcessor postProcessor_;

  //! Free surface integrator module
  solver::FreeSurfaceIntegrator freeSurfaceIntegrator_;

  //! Free surface writer module
  writer::FreeSurfaceWriter freeSurfaceWriter_;

  //! Analysis writer module
  writer::AnalysisWriter analysisWriter_;

  //! Wavefield output module
  writer::WaveFieldWriter waveFieldWriter_;

  //! Fault output module
  writer::FaultWriter faultWriter_;

  //! Receiver writer module
  writer::ReceiverWriter receiverWriter_;

  //! Energy writer module
  writer::EnergyOutput energyOutput_;

  //! Flop Counter
  monitoring::FlopCounter flopCounter_;

  //! TimeMirror Managers
  std::pair<seissol::ITM::InstantaneousTimeMirrorManager,
            seissol::ITM::InstantaneousTimeMirrorManager>
      timeMirrorManagers_;

  //! time stamp which can be used for backuping files of previous runs
  std::string backupTimeStamp_;

  std::optional<std::string> checkpointLoadFile_;

  std::optional<std::size_t> executionPlaceCutoff_;

  utils::Env env_;

  double timestepScale_{1.0};

  public:
  SeisSol(const initializer::parameters::SeisSolParameters& parameters, const utils::Env& env)
      : outputManager_(*this), seissolParameters_(parameters),
        memoryManager_(std::make_unique<initializer::MemoryManager>(*this)), timeManager_(*this),
        freeSurfaceWriter_(*this), analysisWriter_(*this), waveFieldWriter_(*this),
        faultWriter_(*this), receiverWriter_(*this), energyOutput_(*this),
        timeMirrorManagers_(*this, *this), env_(env) {}

  SeisSol(const SeisSol&) = delete;
  SeisSol(SeisSol&&) = delete;
  auto operator=(const SeisSol&) = delete;
  auto operator=(SeisSol&&) = delete;
};

} // namespace seissol

#endif // SEISSOL_SRC_SEISSOL_H_
