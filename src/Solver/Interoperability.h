/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015-2020, SeisSol Group
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

#include <memory>
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>
#include <Initializer/typedefs.hpp>
#include <SourceTerm/NRF.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>
#include <Physics/InitialField.h>
#include "Equations/datastructures.hpp"


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
    TravellingWaveParameters m_travellingWaveParameters;
    
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

    std::vector<Eigen::Vector3d>           m_recPoints;

    //! Vector of initial conditions
    std::vector<std::unique_ptr<physics::InitialField>> m_iniConds;


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

   void setTravellingWaveInformation(const double* origin, const double* kVec, const double* ampField);

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
   void initializeClusteredLts(int clustering, bool enableFreeSurfaceIntegration, bool usePlasticity);
   void initializeMemoryLayout(int clustering, bool enableFreeSurfaceIntegration, bool usePlasticity);
   void bindFaultOutputManager();

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
   //! \todo Documentation
   void setupNRFPointSources( char const* fileName );
#endif

   //! \todo Documentation
   void setupFSRMPointSources( double const* momentTensor,
                               double const* solidVelocityComponent,
                               double const* pressureComponent,
                               double const* fluidVelocityComponent,
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
                          bool    anelasticity,
                          bool    plasticity,
                          bool    anisotropy,
                          bool    poroelasticity,
                          double* materialVal,
                          double* bulkFriction,
                          double* plastCo,
                          double* iniStress,
                          double* waveSpeeds );

    void fitAttenuation(  double rho,
                          double mu,
                          double lambda,
                          double Qp,
                          double Qs,
                          seissol::model::ViscoElasticMaterial& material );

   /**
    * Adds a receiver at the specified location.
    *
    * @param x,y,z coordinates in physical space
    **/
   void addRecPoint(double x, double y, double z) {
     m_recPoints.emplace_back(Eigen::Vector3d(x, y, z));
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

   /***
    * Get the wavespeeds for elastic materials for a given cell
    **/
   void getWaveSpeeds( double* i_materialVal,
                       int     i_numMaterialVals,
                       double* o_waveSpeeds );
   /**
    * Sets the initial loading for a cell (plasticity).
    *
    * @param i_meshId mesh id.
    * @param i_initialLoading initial loading (stress tensor).
    **/
   void setInitialLoading( int    i_meshId,
                           double *i_initialLoading );

   /**
    * Sets the parameters for a cell (plasticity).
    *
    * @param i_meshId mesh id.
    * @param i_plasticParameters cell dependent plastic Parameters (volume, cohesion...).
    **/
   void setPlasticParameters( int    i_meshId,
                              double *i_plasticParameters );

   void setTv(double tv);

   /**
    * \todo Move this somewhere else when we have a C++ main loop.
    **/
   void initializeCellLocalMatrices(bool usePlasticity);

   template<typename T>
   void synchronize(seissol::initializers::Variable<T> const& handle);

   /**
    * Synchronizes the cell local material data.
    **/
   void synchronizeCellLocalData(bool usePlasticity);

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
			int numSides, int numBndGP, int refinement, int* outputMask, int* plasticityMask,
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
   std::vector<std::unique_ptr<physics::InitialField>> const& getInitialConditions() {
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


   /*
    * Not in use any more
    * TODO(Ravil, Sebastian): remove this and all callees inside
    */
   void evaluateFrictionLaw(  int face,
                              real QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()],
                              real QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()],
                              real   imposedStatePlus[seissol::tensor::QInterpolated::size()],
                              real   imposedStateMinus[seissol::tensor::QInterpolated::size()],
                              double i_fullUpdateTime,
                              double timePoints[CONVERGENCE_ORDER],
                              double timeWeights[CONVERGENCE_ORDER],
                              seissol::model::IsotropicWaveSpeeds const& waveSpeedsPlus,
                              seissol::model::IsotropicWaveSpeeds const& waveSpeedsMinus );


  /**
  * Code added by Adrian
   * Get function are only required as long as the dynRup Initializtation is still written in Fortran
  */

  /**
   * get initial values from fortran
   * for each ltsFace mapped to the corresponding fortran mesh face Dynamic Rupture
   * used in ltsFace loop to initialize all missing parameters in initializers::DynamicRupture.h
   *
   * @param ltsFace current ltsFace to get Parameters
   * @param meshFace corresponding meshFace (indexing in fortran) to get Parameters saved in DRFaceInformation[ltsFace].meshFace
   * @param initialStressInFaultCS gets initialStressInFaultCS
   * @param mu gets initial mu
   * @param slipRate1 gets initial sliprate in direction 1
   * @param slipRate2 gets initial sliprate in direction 2
   * @param RF gets intial value (bool) of rupture front output per GP
   **/
  void getDynRupParameters(int ltsFace, unsigned meshFace, real (*initialStressInFaultCS)[init::QInterpolated::Stop[0]][6],
      real (*mu)[seissol::init::QInterpolated::Stop[0]], real  (*slipRate1)[init::QInterpolated::Stop[0]],
      real (*slipRate2)[init::QInterpolated::Stop[0]], bool  (*RF)[ init::QInterpolated::Stop[0] ] );

  /**
   * get initial values from fortran
   * for each ltsFace mapped to the corresponding fortran mesh face Dynamic Rupture
   *
   * @param stateVar    State variable used at Rate-and-state friction laws, gets EQN%IniStateVar
   *
   **/
  void getDynRupStateVar(int ltsFace, unsigned meshFace, real (*stateVar)[init::QInterpolated::Stop[0]]);


  /**
   * get initial values from fortran
   * for each ltsFace mapped to the corresponding fortran mesh face Dynamic Rupture
   *
   * @param nucleationStressInFaultCS gets nucleationStressInFaultCS
   *
   **/
  void getDynRupNucStress(int ltsFace, unsigned meshFace, real (*nucleationStressInFaultCS)[init::QInterpolated::Stop[0]][6]);



  /**
 * get initial values from fortran
 * for each ltsFace mapped to the corresponding fortran mesh face Dynamic Rupture
 * used in ltsFace loop to initialize all missing parameters in initializers::DynamicRupture.h for FL = 2
 *
 * @param ltsFace current ltsFace to get Parameters
 * @param meshFace corresponding meshFace (indexing in fortran) to get Parameters in DRFaceInformation[ltsFace].meshFace
 * @param i_RS_a      RS constitutive parameter "a"
 * @param i_RS_sl0    Reference slip
 * @param i_RS_sr0    Reference slip rate
 **/
  void getDynRupFL_3(int ltsFace,  unsigned meshFace,
                                              real *i_RS_a,
                                              real *i_RS_sl0,
                                              real *i_RS_sr0);

  /**
   * get initial values from fortran
   * for FL103 Thermal Pressure (TP)
   *
   *
   * @param TP_grid     grid for TP
   * @param TP_DFinv    inverse Fourier coefficients
   **/
  void getDynRupTP(real TP_grid[seissol::dr::numberOfTPGridPoints], real TP_DFinv[seissol::dr::numberOfTPGridPoints]);

  void copyFrictionOutputInitialStressInFaultCS(unsigned numberOfCells, real (*initialStressInFaultCS)[init::QInterpolated::Stop[0]][6]);

  /**
   * Temporary Interoperability function for Dynamic rupture outputs
   * copy values from C++ computation back to Fortran output writer.
   * copy parameters, which are present in all friction laws
   **/
  void copyFrictionOutputToFortranGeneral(unsigned ltsFace,
                                          unsigned meshFace,
                                          real  (*slip)[init::QInterpolated::Stop[0]],
                                          real  (*slipStrike)[init::QInterpolated::Stop[0]],
                                          real  (*slipDip)[init::QInterpolated::Stop[0]],
                                          real  (*ruptureTime)[init::QInterpolated::Stop[0]],
                                          real  (*dynStressTime)[init::QInterpolated::Stop[0]],
                                          real  (*peakSlipRate)[init::QInterpolated::Stop[0]],
                                          real  (*tractionXY)[init::QInterpolated::Stop[0]],
                                          real  (*tractionXZ)[init::QInterpolated::Stop[0]]
  );

  /**
   * Temporary Interoperability function for Dynamic rupture outputs
   * copy values from C++ computation back to Fortran output writer.
   * copy parameters, which differ from friction law to friction law
   **/
  void copyFrictionOutputToFortranSpecific(unsigned ltsFace,
                                           unsigned meshFace,
                                           real *averagedSlip,
                                           real (*slipRateStrike)[init::QInterpolated::Stop[0]],
                                           real (*slipRateDip)[init::QInterpolated::Stop[0]],
                                           real (*mu)[seissol::init::QInterpolated::Stop[0]]
  );

  void copyFrictionOutputToFortranStateVar(unsigned ltsFace, unsigned meshFace,
                                      real  (*stateVar)[init::QInterpolated::Stop[0]]
  );

  void copyFrictionOutputToFortranStrength(unsigned ltsFace, unsigned meshFace,
                                           real  (*strength)[init::QInterpolated::Stop[0]]
  );

  void copyFrictionOutputToFortranInitialStressInFaultCS(unsigned ltsFace, unsigned meshFace,
                                                         real  (*initialStressInFaultCS)[init::QInterpolated::Stop[0]][6],
                                                         std::vector<std::array<real, init::QInterpolated::Stop[0]>>& iniBulkXX,
                                                         std::vector<std::array<real, init::QInterpolated::Stop[0]>>& iniBulkYY,
                                                         std::vector<std::array<real, init::QInterpolated::Stop[0]>>& iniBulkZZ,
                                                         std::vector<std::array<real, init::QInterpolated::Stop[0]>>& iniShearXY,
                                                         std::vector<std::array<real, init::QInterpolated::Stop[0]>>& iniShearYZ,
                                                         std::vector<std::array<real, init::QInterpolated::Stop[0]>>& iniShearXZ);

  void initializeFaultOutput();


   /**
    * Prepare element wise faultoutput
    *
    * @param time The current simulation time
    */
   void calcElementwiseFaultoutput( double time );

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

   /**
    * reports memory consumed by each device i.e., GPUs
    */
   void reportDeviceMemoryStatus();

   /**
    * Deallocates memory manager
    */
   void deallocateMemoryManager();
};

#endif
