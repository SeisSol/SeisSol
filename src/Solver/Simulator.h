// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

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
     * update m_abort to abort the main loop of the simulation (see terminator)
     */
    void abort();

    /**
     * Simulates until finished.
     **/
    void simulate(seissol::SeisSol& seissolInstance);
};


#endif // SEISSOL_SRC_SOLVER_SIMULATOR_H_

