// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_SIMULATION_SEISSOL_H_
#define SEISSOL_SRC_SIMULATION_SEISSOL_H_

#include "Alignment.h"
#include "Core/Common/Executor.h"
#include "Core/IO/AsyncIO.h"
#include "Core/IO/Manager.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Modules/Modules.h"
#include "Core/Monitoring/EpochSynchronization.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/MiniSeisSolWriter.h"
#include "Core/Monitoring/RunSummary.h"
#include "Core/Parallel/Pin.h"
#include "Core/Parallel/SystemInfo.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Clustering/VertexWeights/MeasuredCellCosts.h"
#include "Core/TimeStepping/Feature.h"
#include "Core/TimeStepping/TimeManager.h"
#include "Features/InitialConditions/Output/AnalysisWriter.h"
#include "Features/SourceTerm/Setup/Manager.h"
#include "Features/TimeMirror/Solver/InstantaneousTimeMirrorManager.h"
#include "Simulation/EnergyOutput.h"
#include "Simulation/MemoryManager.h"
#include "Simulation/Parameters/SeisSolParameters.h"
#include "Simulation/Simulator.h"
#include "WavePropagation/Output/DofSync.h"
#include "WavePropagation/Output/FreeSurfaceIntegrator.h"
#include "WavePropagation/Output/ReceiverWriter.h"

#include <cstdint>
#include <list>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/env.h>
#include <utils/logger.h>
#include <vector>

namespace seissol {

namespace geometry {
class MeshReader;
} // namespace geometry

/**
 * Subobjects which every member of SeisSol relies on.
 *
 * Base class subobjects are initialized before any non-static data member and destroyed after
 * all of them. Every member of SeisSol therefore finds these set up, while it is constructed as
 * well as while it is destroyed, whatever order the members are declared in.
 *
 * Pinning has to be set up before the first OpenMP call: MPI sets the affinity mask for the
 * process, and the OpenMP runtime replaces it with the pinning given by e.g. OMP_PLACES as soon
 * as it starts up. That is why it comes first.
 *
 * The members read their settings from the environment, some of them already in their
 * constructor.
 *
 * Modules holds non-owning pointers to the modules registered with it. All of those are members
 * of SeisSol.
 */

/**
 * @todo Initialize rank
 */
class SeisSol : public RunContext {
  public:
  /**
   * Cleanup data structures
   */
  ~SeisSol() override;

  /**
   * Initialize C++ part of the program
   */
  bool init();

  /**
   * Finalize SeisSol; and MPI with it, unless another instance runs after this one (see
   * discretization.vertexweightcycles).
   */
  void finalize(bool finalizeMpi = true);

  /**
   * The costs of the cells measured in the run before (in this process), which partition the mesh
   * of this run; none by default. They have to outlive the initialization.
   */
  void setMeasuredCosts(const initializer::MeasuredCellCosts* costs) { measuredCosts_ = costs; }
  [[nodiscard]] const initializer::MeasuredCellCosts* measuredCosts() const {
    return measuredCosts_;
  }

  /**
   * Makes this run a short one that measures what the cells cost, for the partition of the next
   * one: it runs `steps` time steps of the largest time cluster, and keeps the costs (cellCosts()).
   */
  void setMeasuringSteps(unsigned steps) { measuringSteps_ = steps; }
  [[nodiscard]] std::optional<unsigned> measuringSteps() const { return measuringSteps_; }

  //! what the cells cost in this run, where it measured that
  std::optional<initializer::MeasuredCellCosts>& cellCosts() { return cellCosts_; }

  void loadCheckpoint(const std::string& file);

  Executor executionPlace(std::size_t clusterSize) {
    constexpr auto DefaultDevice = isDeviceOn() ? Executor::Device : Executor::Host;
    if (executionPlaceCutoff_.has_value()) {
      // An empty cluster computes nothing anywhere. On the host, it would only make its neighbors
      // exchange data with the host, and keep super-timesteps from getting recorded.
      if (executionPlaceCutoff_.value() <= clusterSize || clusterSize == 0) {
        return DefaultDevice;
      } else {
        return Executor::Host;
      }
    } else {
      return DefaultDevice;
    }
  }

  void setExecutionPlaceCutoff(std::size_t size);

  initializer::MemoryManager& memoryManager() { return *memoryManager_; }

  solver::TimeManager& timeManager() { return timeManager_; }

  Simulator& simulator() { return simulator_; }

  solver::FreeSurfaceIntegrator& freeSurfaceIntegrator() { return freeSurfaceIntegrator_; }

  writer::AnalysisWriter& analysisWriter() { return analysisWriter_; }

  /**
   * Get the receiver writer module
   */
  /// A new writer of an output of receivers, which the run keeps until its end.
  writer::ReceiverWriter& addReceiverWriter() { return receiverWriters_.emplace_back(*this); }

  /**
   * Get the receiver writer module
   */

  /**
   * Get the receiver writer module
   */
  writer::DofSync& dofSync() { return dofSync_; }

  /**
   * Get the energy writer module
   */
  writer::EnergyOutput& energyOutput() { return energyOutput_; }

  /**
   * Get the flop counter
   */

  /**
   * Get the summary of the run, which it writes for programs to read
   */
  monitoring::RunSummary& runSummary() { return runSummary_; }

  /**
   * Get the synchronization points that end the simulation epochs of their own
   */
  monitoring::EpochSynchronization& epochSynchronization() { return epochSynchronization_; }

  /**
   * Get the mini SeisSol writer; it holds the measurement until the output is set up
   */
  writer::MiniSeisSolWriter& miniSeisSolWriter() { return miniSeisSolWriter_; }

  const std::optional<std::string>& checkpointLoadFile() { return checkpointLoadFile_; }
  /**
   * Reference for timeMirrorManagers to be accessed externally when required
   */
  std::list<seissol::physics::InstantaneousTimeMirrorManager>& getTimeMirrorManagers() {
    return timeMirrorManagers_;
  }

  /// The features of the run (see solver::Feature), in the order their hooks run in.
  std::vector<std::unique_ptr<solver::Feature>>& features() { return features_; }

  /**
   * Set the mesh reader
   */
  void setMeshReader(std::unique_ptr<seissol::geometry::MeshReader> meshReader);

  /**
   * Get the mesh reader
   */
  const seissol::geometry::MeshReader& meshReader() const { return *meshReader_; }

  /**
   * Get the mesh reader
   */
  seissol::geometry::MeshReader& meshReader() { return *meshReader_; }

  const seissol::initializer::parameters::SeisSolParameters& parameters() const {
    return seissolParameters_;
  }

  /**
   * Deletes memoryManager. MemoryManager desctructor will destroy all storage structures and
   * memoryAllocator i.e., the main components of SeisSol. Therefore, call this function
   * at the very end of a program execution
   */
  void deleteMemoryManager() { memoryManager_.reset(nullptr); }

  /*
   * sets a time stamp for backuping
   * */
  void setBackupTimeStamp(const std::string& stamp);

  void setTimestepScale(double scale) { timestepScale_ = scale; }

  double timestepScale() const { return timestepScale_; }

  /*
   * The rate vector the clustering actually settled on, normalized to full length.
   *
   * This is not necessarily the vector from the parameter file: a search that chooses the
   * ladder rather than only the wiggle factor publishes its result here, and everything
   * downstream of the mesh reader has to use this one. Set alongside the timestep scale,
   * which is the other half of the same decision.
   * */
  void setEffectiveClusterRates(const std::vector<std::uint64_t>& rates) {
    effectiveClusterRates_ = rates;
  }

  const std::vector<std::uint64_t>& getEffectiveClusterRates() const {
    return effectiveClusterRates_;
  }

  /*
   * returns the backup time stamp
   * */
  const std::string& backupTimeStamp() const { return backupTimeStamp_; }

  /**
   * What the run knows about the nodes and devices it runs on.
   */
  const SystemInfo& systemInfo() const { return systemInfo_; }

  private:
  //! Collection of Parameters
  seissol::initializer::parameters::SeisSolParameters seissolParameters_;

  //! Async I/O handler (needs to be initialize before other I/O modules)

  //! Mesh Reader
  std::unique_ptr<seissol::geometry::MeshReader> meshReader_;

  //! Memory Manager
  std::unique_ptr<initializer::MemoryManager> memoryManager_{nullptr};

  //! Time Manager
  solver::TimeManager timeManager_;

  //! Simulator
  Simulator simulator_;

  //! Free surface integrator module
  solver::FreeSurfaceIntegrator freeSurfaceIntegrator_;

  //! Analysis writer module
  writer::AnalysisWriter analysisWriter_;

  //! Receiver (off-fault) writer module
  std::list<writer::ReceiverWriter> receiverWriters_;

  //! Receiver (on-fault) writer module

  //! Dof sync (for output) module
  writer::DofSync dofSync_;

  //! Energy writer module
  writer::EnergyOutput energyOutput_;

  //! Flop Counter

  //! Run summary
  monitoring::RunSummary runSummary_;

  //! Synchronization points that end the simulation epochs of their own
  monitoring::EpochSynchronization epochSynchronization_;

  //! Mini SeisSol measurement, written with the rest of the output
  writer::MiniSeisSolWriter miniSeisSolWriter_;

  //! TimeMirror Managers
  std::list<seissol::physics::InstantaneousTimeMirrorManager> timeMirrorManagers_;

  //! the features of the run, made with its layout
  std::vector<std::unique_ptr<solver::Feature>> features_;

  //! time stamp which can be used for backuping files of previous runs
  std::string backupTimeStamp_;

  std::optional<std::string> checkpointLoadFile_;

  const initializer::MeasuredCellCosts* measuredCosts_{nullptr};
  std::optional<unsigned> measuringSteps_;
  std::optional<initializer::MeasuredCellCosts> cellCosts_;

  std::optional<std::size_t> executionPlaceCutoff_;

  double timestepScale_{1.0};
  std::vector<std::uint64_t> effectiveClusterRates_;

  const SystemInfo& systemInfo_;

  public:
  SeisSol(initializer::parameters::SeisSolParameters parameters,
          const utils::Env& env,
          const SystemInfo& systemInfo)
      : RunContext(env), seissolParameters_(std::move(parameters)),
        memoryManager_(std::make_unique<initializer::MemoryManager>(*this)),
        timeManager_(
            *this,
            [this](initializer::AllocationPlace place) { memoryManager_->synchronizeTo(place); }),
        analysisWriter_(*this), energyOutput_(*this), systemInfo_(systemInfo) {}

  SeisSol(const SeisSol&) = delete;
  SeisSol(SeisSol&&) = delete;
  auto operator=(const SeisSol&) = delete;
  auto operator=(SeisSol&&) = delete;
};

} // namespace seissol

#endif // SEISSOL_SRC_SIMULATION_SEISSOL_H_
