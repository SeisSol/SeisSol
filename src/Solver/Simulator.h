// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_SOLVER_SIMULATOR_H_
#define SEISSOL_SRC_SOLVER_SIMULATOR_H_

namespace seissol {
  class Simulator;
  class SeisSol;
}

/**
 * Simulator, which takes care of the simulation: Sync. times, wave field output, checkpoints.
 **/
class seissol::Simulator {
  // private:
    //! current time of the simulation
    double m_currentTime;

    //! final time of the simulation
    double m_finalTime;

    //! usePlasticity = true if plasticity is on
    bool m_usePlasticity;

    //! last time a checkpoint was written
    double m_checkPointTime;

    //! time interval of the checkpoints
    double m_checkPointInterval;

    //! If true, a checkpoint is loaded before the simulation
    bool m_loadCheckPoint;

    //! If true, the while loop of the simulation will be aborted (see terminator)
    bool m_abort;

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
     * @param i_finalTime final time.
     **/
    void setFinalTime( double i_finalTime );

    /**
     * Sets the m_usePlasticity
     *
     * @param i_plasticity = 1 if plasticity is on
     **/
    void setUsePlasticity( bool plasticity );

    /**
     * Sets the current time of the simulation (useful for checkpoints)
     *
     * @param i_currentTime current time
     */
    void setCurrentTime( double i_currentTime );

    /**
     * Activates checkpoint loading at the beginning of the simulation
     */
    void loadCheckPoint();

    /**
     * Sets the interval for checkpointing.
     *
     * @param i_checkPointInterval check point interval.
     **/
    void setCheckPointInterval( double i_checkPointInterval );

    /**
     * Returns if the simulator is going to write check points.
     */
    bool checkPointingEnabled();

    /**
     * update m_abort to abort the main loop of the simulation (see terminator)
     */
    void abort();

    /**
     * Simulates until finished.
     **/
    void simulate(seissol::SeisSol& seissolInstance);
};


#endif // SEISSOL_SRC_SOLVER_SIMULATOR_H_

