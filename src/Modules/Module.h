/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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

#ifndef MODULE_H
#define MODULE_H

namespace seissol {

/**
 * Base class for all modules
 */
class Module {
  private:
  /** The synchronization interval for this module */
  double isyncInterval;

  /** The next synchronization point for this module */
  double nextSyncPoint;

  /** The last time when syncPoint was called */
  double lastSyncPoint;

  public:
  Module();

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
  double syncInterval() const { return isyncInterval; }

  /**
   * Set the synchronization interval for this module
   *
   * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
   */
  void setSyncInterval(double interval);
};

} // namespace seissol

#endif // MODULE_H
