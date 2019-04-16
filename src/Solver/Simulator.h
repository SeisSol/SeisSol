/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2017, SeisSol Group
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
 * Entry point of the simulation.
 **/

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

namespace seissol {
  class Simulator;
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

    //! last time a checkpoint was written
    double m_checkPointTime;

    //! time interval of the checkpoints
    double m_checkPointInterval;

    //! last time energies were printed
    double m_printEnergiesTime;

    //! time interval when to print energies
    double m_printEnergiesInterval;

    //! If true, a checkpoint is loaded before the simulation
    bool m_loadCheckPoint;

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
     * Sets the interval for printing Energies.
     **/
    void setPrintEnergiesInterval( double i_printEnergiesInterval );

    /**
     * Returns if the simulator is going to write check points.
     */
    bool checkPointingEnabled();

    /**
     * Simulates until finished.
     **/
    void simulate();
};

#endif
