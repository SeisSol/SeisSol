// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_MODULES_MODULE_H_
#define SEISSOL_SRC_MODULES_MODULE_H_

namespace seissol {

/**
 * Base class for all modules
 */
class Module {
  private:
  /** The synchronization interval for this module */
  double isyncInterval{0};

  /** The next synchronization point for this module */
  double nextSyncPoint{0};

  /** The last time when syncPoint was called */
  double lastSyncPoint;

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
  virtual void simulationStart() {}

  virtual void simulationEnd() {}

  virtual void shutdown() {}

  /**
   * Called at synchronization points
   */
  virtual void syncPoint(double currentTime) {}

  protected:
  [[nodiscard]] double syncInterval() const { return isyncInterval; }

  /**
   * Set the synchronization interval for this module
   *
   * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
   */
  void setSyncInterval(double interval);
};

} // namespace seissol

#endif // SEISSOL_SRC_MODULES_MODULE_H_
