// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MODULES_MODULE_H_
#define SEISSOL_SRC_MODULES_MODULE_H_

#include "Memory/HostData.h"

#include <optional>
namespace seissol {

/**
 * Base class for all modules
 */
class Module {
  private:
  /** The synchronization interval for this module */
  double isyncInterval_{0};

  /** The next synchronization point for this module */
  double nextSyncPoint_{0};

  /** The last time when syncPoint was called */
  double lastSyncPoint_;

  public:
  Module();

  virtual ~Module();

  /**
   * Called by {@link Modules} at every synchronization point
   *
   * We have to ensure that this is "our" synchronization point before
   * calling {@link syncPoint}.
   *
   * @return The next synchronization point for this module
   */
  double potentialSyncPoint(double currentTime, double timeTolerance, bool forceSyncPoint);

  /**
   * Whether potentialSyncPoint() would call {@link syncPoint} at `currentTime`; without side
   * effects.
   */
  [[nodiscard]] bool dueAt(double currentTime, double timeTolerance, bool forceSyncPoint) const;

  /**
   * Notes that `currentTime` was a synchronization point of this module, as potentialSyncPoint()
   * does after calling {@link syncPoint}.
   */
  void markSyncPoint(double currentTime);

  /** The next synchronization point of this module. */
  [[nodiscard]] double nextSyncPoint() const { return nextSyncPoint_; }

  /**
   * Whether the second part of a synchronization point (processSyncPoint()) may run on the thread
   * of the outputs (see Modules), while the device goes on with the next epoch. Then
   * captureSyncPoint() takes all that changes in that epoch, and processSyncPoint() reads only
   * what it took, and the host copies of the data (see hostData()).
   */
  [[nodiscard]] virtual bool processesAsynchronously() const { return false; }

  /**
   * The variables of the storages this module reads on the host at its synchronization points,
   * which get copied there from the device before; all of them, unless the module says otherwise.
   */
  [[nodiscard]] virtual initializer::HostData hostData() const {
    return initializer::HostData::everything();
  }

  /**
   * Called by {@link Modules} before the simulation starts to set the synchronization point.
   *
   * This is only called for modules that register for the SYNCHRONIZATION_POINT hook.
   */
  void setSimulationStartTime(double time);

  //
  // Potential hooks
  //

  /**
   * Called before initializing MPI
   */
  virtual void preMPI() {}

  /**
   * Called after MPI initialization
   */
  virtual void postMPIInit() {}

  /**
   * Called before mesh initialization
   */
  virtual void preMesh() {}

  /**
   * Called after mesh initialization
   */
  virtual void postMesh() {}

  /**
   * Called before LTS initialization
   */
  virtual void preLtsInit() {}

  /**
   * Called after LTS initialization
   */
  virtual void postLtsInit() {}

  /**
   * Called before the model is initialized
   */
  virtual void preModel() {}

  /**
   * Called after the model is initialized
   */
  virtual void postModel() {}

  /**
   * Called before the actual simulation.
   *
   * Only called when the simulation is not started from a checkpoint
   */
  virtual void simulationStart(std::optional<double> checkpointTime) {}

  virtual void simulationEnd() {}

  virtual void shutdown() {}

  /**
   * Called at synchronization points; by default, the two parts of one, on the calling thread.
   */
  virtual void syncPoint(double currentTime) {
    captureSyncPoint(currentTime);
    processSyncPoint(currentTime);
  }

  /**
   * The first part of a synchronization point, which runs on the main thread, with the data just
   * copied to the host, before the next epoch starts.
   */
  virtual void captureSyncPoint(double /*currentTime*/) {}

  /**
   * The second part of a synchronization point; on the thread of the outputs if
   * processesAsynchronously().
   */
  virtual void processSyncPoint(double /*currentTime*/) {}

  protected:
  [[nodiscard]] double syncInterval() const { return isyncInterval_; }

  /**
   * Set the synchronization interval for this module
   *
   * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
   */
  void setSyncInterval(double interval);
};

} // namespace seissol

#endif // SEISSOL_SRC_MODULES_MODULE_H_
