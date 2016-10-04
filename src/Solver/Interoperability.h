/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 * C++/Fortran-interoperability.
 **/

#ifndef INTEROPERABILITY_H_
#define INTEROPERABILITY_H_

#include <vector>
#include <Initializer/typedefs.hpp>
#include <Kernels/Time.h>
#include <SourceTerm/NRF.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>

namespace seissol {
  class Interoperability;
}

/**
 *  C++ binding to all required fortran functions.
 **/
class seissol::Interoperability {
  // private:
    //! time kernel
    seissol::kernels::Time m_timeKernel;

    /* Brain dump of SeisSol's Fortran parts:
     * Raw fotran-pointer to cope with limited modularity of the
     * source, receiver and dynamic rupture functions.
     */
    //! domain (tDomain)
    void *m_domain;

    /*
     * Brain dump of SeisSol's C parts.
     */
    //! cluster-local mesh structure
    struct MeshStructure *m_meshStructure;

    //! time stepping
    struct TimeStepping m_timeStepping;

    //! global data
    struct GlobalData* m_globalData;

    seissol::initializers::LTSTree*   m_ltsTree;
    seissol::initializers::LTS*       m_lts;

    //! Lookup table relating mesh to cells
    seissol::initializers::Lut        m_ltsLut;

 public:
   /**
    * Constructor.
    **/
   Interoperability();

   ~Interoperability();

   /**
    * Sets the fortran domain.
    *
    * @param i_domain domain.
    */
   void setDomain( void *i_domain );

   /**
    * Sets the time step width of a cell.
    *
    * @param i_meshId mesh id.
    * @param i_timeStepWidth time step width of the cell, Fortran notation is assumed (starting at 0).
    **/
   void setTimeStepWidth( int    i_meshId,
                          double i_timeStepWidth );

   /**
    * Initializes clustered local time stepping.
    *   1) Derivation of the LTS layout
    *   2) Memory Setup
    *   3) Generation of LTS computer clusters
    *
    * Clustering stategy is mapped as follows:
    *   1:  Global time stepping
    *   2+: Fixed rate between clusters
    *
    * @param i_clustering clustering strategy
    **/
   void initializeClusteredLts( int i_clustering );

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
   //! \todo Documentation
   void setupNRFPointSources( char const* fileName );
#endif

   //! \todo Documentation
   void setupFSRMPointSources( double const* momentTensor,
                               int           numberOfSources,
                               double const* centres,
                               double const* strikes,
                               double const* dips,
                               double const* rakes,
                               double const* onsets,
                               double const* areas,
                               double        timestep,
                               int           numberOfSamples,
                               double const* timeHistories );

   /**
    * Adds a receiver at the specified mesh id.
    *
    * @param i_receiverId pointer to the global id of the receiver.
    * @param i_meshId pointer to the mesh id.
    **/
   void addReceiver( int i_receiverId,
                     int i_meshId );

   /**
    * Sets the sampling of the receivers.
    *
    * @param i_receiverSampling sampling of the receivers.
    **/
   void setReceiverSampling( double i_receiverSampling );

   /**
    * Enables dynamic rupture.
    **/
   void enableDynamicRupture();

   /**
    * Set material parameters for cell
    **/
   void setMaterial( int    i_meshId,
                     int    i_side,
                     double* i_materialVal,
                     int    i_numMaterialVals );

   /**
    * Sets the initial loading for a cell (plasticity).
    *
    * @param i_meshId mesh id.
    * @param i_initialLoading initial loading (stress tensor).
    **/
#ifdef USE_PLASTICITY
   void setInitialLoading( int    *i_meshId,
                           double *i_initialLoading );
#endif

   /**
    * Sets the parameters for a cell (plasticity).
    *
    * @param i_meshId mesh id.
    * @param i_plasticParameters cell dependent plastic Parameters (volume, cohesion...).
    **/
#ifdef USE_PLASTICITY
   void setPlasticParameters( int    *i_meshId,
                              double *i_plasticParameters );
#endif

   /**
    * \todo Move this somewhere else when we have a C++ main loop.
    **/
   void initializeCellLocalMatrices();

   template<typename T>
   void synchronize(seissol::initializers::Variable<T> const& handle);

   /**
    * Synchronizes the cell local material data.
    **/
   void synchronizeCellLocalData();

   /**
    * Synchronizes the DOFs in the copy layer.
    **/
   void synchronizeCopyLayerDofs();

   /**
    * Enable wave field plotting.
    *
    * @param i_waveFieldInterval plotting interval of the wave field.
    * @param i_waveFieldFilename file name prefix of the wave field.
    **/
   void enableWaveFieldOutput( double i_waveFieldInterval, const char *i_waveFieldFilename );

   /**
    * Enable checkpointing.
    *
    * @param i_checkPointInterval check pointing interval.
    * @param i_checkPointFilename file name prefix for checkpointing.
    * @param i_checkPointBackend name of the checkpoint backend
    **/
   void enableCheckPointing( double i_checkPointInterval,
		   const char *i_checkPointFilename, const char* i_checkPointBackend );

   void initializeIO(double* mu, double* slipRate1, double* slipRate2,
			  double* slip, double* slip1, double* slip2, double* state, double* strength,
			  int numSides, int numBndGP, int refinement, int* outputMask, double* outputRegionBounds);

   /**
    * Get the current dynamic rupture time step
    *
    * @param o_timeStep The dynamic rupture time step
    */
   void getDynamicRuptureTimeStep(int &o_timeStep);

   /**
    * Adds the specified update to dofs.
    *
    * @param i_mesh mesh id of the cell, Fortran notation is assumed - starting at 1 instead of 0.
    * @param i_update update which is applied.
    **/
   void addToDofs( int      i_meshId,
                   double*  i_update,
                   int      numberOfQuantities );

   /**
    * Writes the receivers.
    *
    * @param i_fullUpdateTime full update time of the DOFs relevant to the receivers.
    * @param i_timeStepWidth time step width of the next update.
    * @param i_receiverTime time of the receivers / last write.
    * @param i_receiverIds global ids of the receivers.
    **/
   void writeReceivers( double              i_fullUpdateTime,
                        double              i_timeStepWidth,
                        double              i_receiverTime,
                        std::vector< int > &i_receiverIds );

   /**
    * Gets the time derivatives (recomputed from DOFs).
    *
    * @param i_meshId mesh id.
    * @param o_timeDerivatives time derivatives in deprecated full storage scheme (including zero blocks).
    **/
   void getTimeDerivatives( int    i_meshId,
                            double  o_timeDerivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] );

   /**
    * Gets the time derivatives and integrated DOFs of two face neighbors.
    *
    * @param i_meshId mesh id.
    * @param i_localFaceId local id of the face neighbor.
    * @param i_timeStepWidth time step width used for the time integration.
    * @param o_timeDerivativesCell time derivatives of the cell in deprecated full storage scheme (including zero blocks).
    * @param o_timeDerivativesNeighbor time derivatives of the cell neighbor in deprecated full storage scheme (including zero blocks).
    * @param o_timeIntegratedCell time integrated degrees of freem of the cell.
    * @param o_timeIntegratedNeighbor time integrated degrees of free of the neighboring cell.
    **/
   void getFaceDerInt( int    _meshId,
                       int    i_localFaceId,
                       double i_timeStepWidth,
                       double o_timeDerivativesCell[CONVERGENCE_ORDER][NUMBER_OF_DOFS],
                       double o_timeDerivativesNeighbor[CONVERGENCE_ORDER][NUMBER_OF_DOFS],
                       double o_timeIntegratedCell[NUMBER_OF_DOFS],
                       double o_timeIntegratedNeighbor[NUMBER_OF_DOFS] );

   /**
    * Gets the DOFs.
    *
    * @param i_meshId mesh id.
    * @param o_dofs degrees of freedom.
    **/
   void getDofs( int    i_meshId,
                 double o_dofs[NUMBER_OF_DOFS] );

   /**
    * Gets the DOFs from the derivatives.
    * Assumes valid storage of time derivatives.
    *
    * @param i_meshId mesh id.
    * @param o_dofs degrees of freedom.
    **/
   void getDofsFromDerivatives( int    i_meshId,
                                double o_dofs[NUMBER_OF_DOFS] );

   /**
    * Gets the neighboring DOFs from the derivatives.
    * Assumes valid storage of time derivatives.
    *
    * @param i_meshId mesh id.
    * @param i_localFaceId local id of the face neighbor.
    * @param o_dofs degrees of freedom.
    **/
   void getNeighborDofsFromDerivatives( int    i_meshId,
                                        int    i_localFaceId,
                                        double o_dofs[NUMBER_OF_DOFS] );

   /**
    * Computes dynamic rupture on the faces.
    *
    * @param i_fullUpdateTime full update time of the respective DOFs.
    * @param i_timeStepWidth time step width of the next full update.
    **/
   void computeDynamicRupture( double i_fullUpdateTime,
                               double i_timeStepWidth );

   /**
    * Computes plasticity.
    *
    * @param i_timeStep time step of the previous update.
    * @param i_plasticParameters cell dependent plasticity parameters
    * @param i_initialLoading initial loading of the associated cell.
    * @param io_dofs degrees of freedom (including alignment).
    * @param io_pstrain plastic strain tensor
    **/
#ifdef USE_PLASTICITY
   void computePlasticity( double   i_timeStep,
		                   double  *i_plasticParameters,
                           double (*i_initialLoading)[NUMBER_OF_BASIS_FUNCTIONS],
                           double  *io_dofs,
						   double  *io_Energy,
						   double  *io_pstrain );
#endif

   /**
    * Computes mInvJInvPhisAtSources[i] = |J|^-1 * M_ii^-1 * phi_i(xi, eta, zeta),
    * where xi, eta, zeta is the point in the reference tetrahedron corresponding to x, y, z.
    *
    * @param x x coordinate
    * @param y y coordinate
    * @param z z coordinate
    * @param element Number of element in that x, y, z resides
    * @param mInvJInvPhisAtSources contains the output
    */
   void computeMInvJInvPhisAtSources( double x,
                                      double y,
                                      double z,
                                      unsigned element,
                                      real mInvJInvPhisAtSources[NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] );

   /**
    * Simulates until the final time is reached.
    *
    * @param i_finalTime final time to reach.
    **/
   void simulate( double i_finalTime );

   /**
    * Finalizes I/O
    */
   void finalizeIO();
};

#endif
