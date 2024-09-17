// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 **/

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

namespace seissol {
class Simulator;
class SeisSol;
} // namespace seissol

/**
 * Simulator, which takes care of the simulation: Sync. times, wave field output, checkpoints.
 **/
class seissol::Simulator {
  // private:
  //! current time of the simulation
  double currentTime;

  //! final time of the simulation
  double finalTime;

  //! usePlasticity = true if plasticity is on
  bool usePlasticity;

  //! last time a checkpoint was written
  double checkPointTime;

  //! time interval of the checkpoints
  double checkPointInterval;

  //! If true, the while loop of the simulation will be aborted (see terminator)
  bool aborted;

  public:
  /**
   * Constructor, which initializes all values.
   * Default:
   *  All times are set to zero.
   *  Wave field and checkpoint are disabled, intervals set to infinity.
   **/
  Simulator();

  /**
   * Sets the final time of the simulation.
   *
   * @param finalTime final time.
   **/
  void setFinalTime(double finalTime);

  /**
   * Sets the usePlasticity
   *
   * @param plasticity = 1 if plasticity is on
   **/
  void setUsePlasticity(bool plasticity);

  /**
   * Sets the current time of the simulation (useful for checkpoints)
   *
   * @param currentTime current time
   */
  void setCurrentTime(double currentTime);

  /**
   * Sets the interval for checkpointing.
   *
   * @param checkPointInterval check point interval.
   **/
  void setCheckPointInterval(double checkPointInterval);

  /**
   * Returns if the simulator is going to write check points.
   */
  bool checkPointingEnabled();

  /**
   * update abort to abort the main loop of the simulation (see terminator)
   */
  void abort();

  /**
   * Simulates until finished.
   **/
  void simulate(seissol::SeisSol& seissolInstance);
};

#endif
