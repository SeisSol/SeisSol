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
 * C++/Fortran-interoperability.
 **/

#ifndef INTEROPERABILITY_H_
#define INTEROPERABILITY_H_

#include <unordered_map>
#include <vector>
#include <glm/vec3.hpp>
#include <Initializer/typedefs.hpp>
#include <SourceTerm/NRF.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>
#include <Physics/InitialField.h>

namespace seissol {
  class Interoperability;
}

/**
 *  C++ binding to all required fortran functions.
 **/
class seissol::Interoperability {
  // private:
    // Type of the initial condition.
    std::string m_initialConditionType;
    
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

    //! Lookup table relating faces to layers
    unsigned*                         m_ltsFaceToMeshFace;
    
    //! Set of parameters that have to be initialized for dynamic rupture
    std::unordered_map<std::string, double*> m_faultParameters;

    std::vector<glm::dvec3>           m_recPoints;

    //! Vector of initial conditions
    std::vector<physics::InitialField*> m_iniConds;

    void initInitialConditions();
 public:
   /**
    * Constructor.
    **/
   Interoperability();

   ~Interoperability();

   /**
    * Sets the type of the initial conditions.
    *
    * @param type The name of the type of the initial conditions.
    */
   void setInitialConditionType(char const *type);

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
    * @param enableFreeSurfaceIntegration
    **/
   void initializeClusteredLts( int i_clustering, bool enableFreeSurfaceIntegration );

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
   //! \todo Documentation
   void setupNRFPointSources( char const* fileName );
#endif

   //! \todo Documentation
   void setupFSRMPointSources( double const* momentTensor,
                               double const* velocityComponent,
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

    //! \todo Documentation
    void initializeModel( char*   materialFileName,
                          int     anelasticity,
                          int     plasticity,
                          double* materialVal,
                          double* bulkFriction,
                          double* plastCo,
                          double* iniStress );

    void fitAttenuation(  double rho,
                          double mu,
                          double lambda,
                          double Qp,
                          double Qs,
                          seissol::model::Material& material );

    void addFaultParameter( std::string const& name,
                           double* memory) {
      m_faultParameters[name] = memory;
    }
    
    //! \todo Documentation
    void initializeFault( char*   modelFileName,
                          int     gpwise,
                          double* bndPoints,
                          int     numberOfBndPoints );

   /**
    * Adds a receiver at the specified location.
    *
    * @param x,y,z coordinates in physical space
    **/
   void addRecPoint(double x, double y, double z) {
     m_recPoints.emplace_back(glm::dvec3(x, y, z));
   }

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
   void setInitialLoading( int    i_meshId,
                           double *i_initialLoading );
#endif

   /**
    * Sets the parameters for a cell (plasticity).
    *
    * @param i_meshId mesh id.
    * @param i_plasticParameters cell dependent plastic Parameters (volume, cohesion...).
    **/
#ifdef USE_PLASTICITY
   void setPlasticParameters( int    i_meshId,
                              double *i_plasticParameters );

   void setTv(double tv);
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
    * Enable free surface output
    */
   void enableFreeSurfaceOutput(int maxRefinementDepth);

   /**
    * Enable checkpointing.
    *
    * @param i_checkPointInterval check pointing interval.
    * @param i_checkPointFilename file name prefix for checkpointing.
    * @param i_checkPointBackend name of the checkpoint backend
    **/
   void enableCheckPointing( double i_checkPointInterval,
		   const char *i_checkPointFilename, const char* i_checkPointBackend );

   /**
    * Gets the integration mask from the parameter file
    **/
   void getIntegrationMask( int* i_integrationMask );

   void initializeIO(double* mu, double* slipRate1, double* slipRate2,
			double* slip, double* slip1, double* slip2, double* state, double* strength,
			int numSides, int numBndGP, int refinement, int* outputMask,
			double* outputRegionBounds,
			double freeSurfaceInterval, const char* freeSurfaceFilename,
      char const* xdmfWriterBackend,
      double receiverSamplingInterval, double receiverSyncInterval);

   /**
    * Copy dynamic rupture variables for output.
    **/
   void copyDynamicRuptureState();

  /**
   * Returns (possibly multiple) initial conditions
   */
   std::vector<physics::InitialField*> const& getInitialConditions() {
     return m_iniConds;
   }

   /**
    * Project initial field on degrees of freedom.
    **/
   void projectInitialField();

   /**
    * Gets the DOFs.
    *
    * @param i_meshId mesh id.
    * @param o_dofs degrees of freedom.
    **/
   void getDofs( int    i_meshId,
                 double o_dofs[tensor::QFortran::size()] );

   /**
    * Gets the DOFs from the derivatives.
    * Assumes valid storage of time derivatives.
    *
    * @param i_meshId mesh id.
    * @param o_dofs degrees of freedom.
    **/
   void getDofsFromDerivatives( int    i_meshId,
                                double o_dofs[tensor::QFortran::size()] );

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
                                        double o_dofs[tensor::QFortran::size()] );
   /**
    * Gets the LTS lookup table.
    */
   seissol::initializers::Lut* getLtsLut();

   /**
    * Gets the type of the initial conditions.
    */
   std::string getInitialConditionType();

   /**
    * Compute fault output.
    *
    * @param i_fullUpdateTime full update time of the respective DOFs.
    * @param i_timeStepWidth time step width of the next full update.
    **/
   void faultOutput( double i_fullUpdateTime, double i_timeStepWidth );

   void evaluateFrictionLaw(  int face,
                              real   godunov[CONVERGENCE_ORDER][seissol::tensor::godunovState::size()],
                              real   imposedStatePlus[seissol::tensor::godunovState::size()],
                              real   imposedStateMinus[seissol::tensor::godunovState::size()],
                              double i_fullUpdateTime,
                              double timePoints[CONVERGENCE_ORDER],
                              double timeWeights[CONVERGENCE_ORDER],
                              seissol::model::IsotropicWaveSpeeds const& waveSpeedsPlus,
                              seissol::model::IsotropicWaveSpeeds const& waveSpeedsMinus );


   /**
    * Prepare element wise faultoutput
    *
    * @param time The current simulation time
    */
   void calcElementwiseFaultoutput( double time );

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
                                      real mInvJInvPhisAtSources[tensor::mInvJInvPhisAtSources::size()] );

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
