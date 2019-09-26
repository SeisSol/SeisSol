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

#include <cstddef>
#include <cstring>

#include "Interoperability.h"
#include "time_stepping/TimeManager.h"
#include "SeisSol.h"
#include <Initializer/CellLocalMatrices.h>
#include <Initializer/InitialFieldProjection.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/time_stepping/common.hpp>
#include <Model/Setup.h>
#include <Monitoring/FlopCounter.hpp>
#include <ResultWriter/common.hpp>

seissol::Interoperability e_interoperability;

/*
 * C bindings
 */
extern "C" {
  // fortran to c
  void c_interoperability_setDomain( void *i_domain ) {
    e_interoperability.setDomain( i_domain );
  }

  void c_interoperability_setTimeStepWidth( int    i_meshId,
                                            double i_timeStepWidth ) {
    e_interoperability.setTimeStepWidth( i_meshId,
                                         i_timeStepWidth );
  }

  void c_interoperability_initializeClusteredLts( int i_clustering, bool enableFreeSurfaceIntegration ) {
    e_interoperability.initializeClusteredLts( i_clustering, enableFreeSurfaceIntegration );
  }

  void c_interoperability_setInitialConditionType(char* type)
  {
    e_interoperability.setInitialConditionType(type);
  }
  
  void c_interoperability_setupNRFPointSources(char* nrfFileName)
  {
#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
    e_interoperability.setupNRFPointSources(nrfFileName);
#endif
  }

  void c_interoperability_setupFSRMPointSources( double*  momentTensor,
                                                 double*  velocityComponent,
                                                 int      numberOfSources,
                                                 double*  centres,
                                                 double*  strikes,
                                                 double*  dips,
                                                 double*  rakes,
                                                 double*  onsets,
                                                 double*  areas,
                                                 double   timestep,
                                                 int      numberOfSamples,
                                                 double*  timeHistories )
  {
    e_interoperability.setupFSRMPointSources( momentTensor,
                                              velocityComponent,
                                              numberOfSources,
                                              centres,
                                              strikes,
                                              dips,
                                              rakes,
                                              onsets,
                                              areas,
                                              timestep,
                                              numberOfSamples,
                                              timeHistories );
  }

  void c_interoperability_initializeModel(  char*   materialFileName,
                                            int     anelasticity,
                                            int     plasticity,
                                            double* materialVal,
                                            double* bulkFriction,
                                            double* plastCo,
                                            double* iniStress )
  {
    e_interoperability.initializeModel( materialFileName,
                                        anelasticity,
                                        plasticity,
                                        materialVal,
                                        bulkFriction,
                                        plastCo,
                                        iniStress );
  }
  
  void c_interoperability_addFaultParameter(  char* name,
                                              double* memory  ) {
    e_interoperability.addFaultParameter(name, memory);
  }

  bool c_interoperability_faultParameterizedByTraction( char* modelFileName ) {
    return seissol::initializers::ParameterDB::faultParameterizedByTraction( std::string(modelFileName) );
  }

  void c_interoperability_initializeFault(  char*   modelFileName,
                                            int     gpwise,
                                            double* bndPoints,
                                            int     numberOfBndPoints ) {    
    e_interoperability.initializeFault(modelFileName, gpwise, bndPoints, numberOfBndPoints);
  }

  void c_interoperability_addRecPoint(double x, double y, double z) {
    e_interoperability.addRecPoint(x, y, z);
  }

  void c_interoperability_enableDynamicRupture() {
    e_interoperability.enableDynamicRupture();
  }

  void c_interoperability_setMaterial( int    i_meshId,
                                       int    i_side,
                                       double* i_materialVal,
                                       int    i_numMaterialVals ) {
    e_interoperability.setMaterial(i_meshId, i_side, i_materialVal, i_numMaterialVals);
  }

#ifdef USE_PLASTICITY
 void c_interoperability_setInitialLoading( int    i_meshId,
                                            double *i_initialLoading ) {
    e_interoperability.setInitialLoading( i_meshId, i_initialLoading );
  }

 void c_interoperability_setPlasticParameters( int    i_meshId,
                                               double *i_plasticParameters ) {
    e_interoperability.setPlasticParameters( i_meshId, i_plasticParameters );
  }

  void c_interoperability_setTv(double tv) {
    e_interoperability.setTv(tv);
  }
#endif

  void c_interoperability_initializeCellLocalMatrices() {
    e_interoperability.initializeCellLocalMatrices();
  }

  void c_interoperability_synchronizeCellLocalData() {
    e_interoperability.synchronizeCellLocalData();
  }

  void c_interoperability_synchronizeCopyLayerDofs() {
    e_interoperability.synchronizeCopyLayerDofs();
  }

  void c_interoperability_enableWaveFieldOutput( double i_waveFieldInterval, const char* i_waveFieldFilename ) {
    e_interoperability.enableWaveFieldOutput( i_waveFieldInterval, i_waveFieldFilename );
  }

  void c_interoperability_enableFreeSurfaceOutput(int maxRefinementDepth) {
	e_interoperability.enableFreeSurfaceOutput(maxRefinementDepth);
  }

  void c_interoperability_enableCheckPointing( double i_checkPointInterval,
		  const char* i_checkPointFilename, const char* i_checkPointBackend ) {
    e_interoperability.enableCheckPointing( i_checkPointInterval,
    		i_checkPointFilename, i_checkPointBackend );
  }

  void c_interoperability_getIntegrationMask( int* i_integrationMask ) {
    e_interoperability.getIntegrationMask( i_integrationMask );
  }

  void c_interoperability_initializeIO( double* mu, double* slipRate1, double* slipRate2,
		  double* slip, double* slip1, double* slip2, double* state, double* strength,
		  int numSides, int numBndGP, int refinement, int* outputMask, double* outputRegionBounds,
		  double freeSurfaceInterval, const char* freeSurfaceFilename, char const* xdmfWriterBackend,
      double receiverSamplingInterval, double receiverSyncInterval) {
	  e_interoperability.initializeIO(mu, slipRate1, slipRate2, slip, slip1, slip2, state, strength,
			numSides, numBndGP, refinement, outputMask, outputRegionBounds,
			freeSurfaceInterval, freeSurfaceFilename, xdmfWriterBackend,
      receiverSamplingInterval, receiverSyncInterval);
  }

  void c_interoperability_projectInitialField() {
    e_interoperability.projectInitialField();
  }

  void c_interoperability_getDofs( int    i_meshId,
                                   double o_timeDerivatives[seissol::tensor::QFortran::size()] ) {
    e_interoperability.getDofs( i_meshId, o_timeDerivatives );
  }

  void c_interoperability_getDofsFromDerivatives( int    i_meshId,
                                                  double o_dofs[seissol::tensor::QFortran::size()] ) {
    e_interoperability.getDofsFromDerivatives( i_meshId, o_dofs );
  }

  void c_interoperability_getNeighborDofsFromDerivatives( int    i_meshId,
                                                          int    i_localFaceId,
                                                          double o_dofs[seissol::tensor::QFortran::size()] ) {
    e_interoperability.getNeighborDofsFromDerivatives( i_meshId, i_localFaceId, o_dofs );
  }

  void c_interoperability_simulate( double i_finalTime ) {
    e_interoperability.simulate( i_finalTime );
  }

  void c_interoperability_finalizeIO() {
	  e_interoperability.finalizeIO();
  }

  // c to fortran
  extern void f_interoperability_computeSource(  void   *i_domain,
                                                 int    *i_meshId,
                                                 double *i_fullUpdateTime,
                                                 double *i_timeStepWidth );

  extern void f_interoperability_copyDynamicRuptureState(void* domain);

  extern void f_interoperability_faultOutput( void   *i_domain,
                                              double *i_fullUpdateTime,
                                              double *i_timeStepWidth );

  extern void f_interoperability_evaluateFrictionLaw( void*   i_domain,
                                                      int     i_face,
                                                      real*   i_godunov,
                                                      real*   i_imposedStatePlus,
                                                      real*   i_imposedStateMinus,
                                                      int     i_numberOfBasisFunctions2D,
                                                      int     i_godunovLd,
                                                      double* i_fullUpdateTime,
                                                      double* timePoints,
                                                      double* timeWeights,
                                                      double  densityPlus,
                                                      double  pWaveVelocityPlus,
                                                      double  sWaveVelocityPlus,
                                                      double  densityMinus,
                                                      double  pWaveVelocityMinus,
                                                      double  sWaveVelocityMinus );

  extern void f_interoperability_calcElementwiseFaultoutput( void *domain,
	                                                     double time );

  extern void f_interoperability_computePlasticity( void    *i_domain,
                                                    double  *i_timestep,
													int    numberOfAlignedBasisFunctions,
													double  *i_plasticParameters,
                                                    double (*i_initialLoading)[NUMBER_OF_BASIS_FUNCTIONS],
                                                    double  *io_dofs,
													double  *io_Energy,
													double  *io_pstrain );

  extern void f_interoperability_computeMInvJInvPhisAtSources( void*    i_domain,
                                                               double   i_x,
                                                               double   i_y,
                                                               double   i_z,
                                                               int      i_elem,
                                                               double*  o_mInvJInvPhisAtSources );

  extern void f_interoperability_fitAttenuation(  void*  i_domain,
                                                  double  rho,
                                                  double  mu,
                                                  double  lambda,
                                                  double  Qp,
                                                  double  Qs,
                                                  double* material );
}

/*
 * C++ functions
 */
seissol::Interoperability::Interoperability() :
  m_initialConditionType(),  m_domain(nullptr), m_ltsTree(nullptr), m_lts(nullptr), m_ltsFaceToMeshFace(nullptr) // reset domain pointer
{
}

seissol::Interoperability::~Interoperability()
{
  delete[] m_ltsFaceToMeshFace;
  for (auto& iniCond : m_iniConds) {
    delete iniCond;
  }
}

void seissol::Interoperability::setInitialConditionType(char const* type) {
  assert(type != nullptr);
  // Note: Pointer to type gets deleted before doing the error computation.
  m_initialConditionType = std::string(type);
}

void seissol::Interoperability::setDomain( void* i_domain ) {
  assert( i_domain != NULL );
  m_domain = i_domain;
}

void seissol::Interoperability::setTimeStepWidth( int    i_meshId,
                                                  double i_timeStepWidth ) {
  seissol::SeisSol::main.getLtsLayout().setTimeStepWidth( (i_meshId)-1, i_timeStepWidth );
}

void seissol::Interoperability::initializeClusteredLts( int i_clustering, bool enableFreeSurfaceIntegration ) {
  // assert a valid clustering
  assert( i_clustering > 0 );

  // either derive a GTS or LTS layout
  if( i_clustering == 1 ) {
    seissol::SeisSol::main.getLtsLayout().deriveLayout( single, 1);
  }
  else {
    seissol::SeisSol::main.getLtsLayout().deriveLayout( multiRate, i_clustering );
  }

  // get the mesh structure
  seissol::SeisSol::main.getLtsLayout().getMeshStructure( m_meshStructure );

  // get time stepping
  seissol::SeisSol::main.getLtsLayout().getCrossClusterTimeStepping( m_timeStepping );


  unsigned* numberOfDRCopyFaces;
  unsigned* numberOfDRInteriorFaces;
  // get cell information & mappings
  seissol::SeisSol::main.getLtsLayout().getDynamicRuptureInformation( m_ltsFaceToMeshFace,
                                                                      numberOfDRCopyFaces,
                                                                      numberOfDRInteriorFaces );

  seissol::SeisSol::main.getMemoryManager().fixateLtsTree(  m_timeStepping,
                                                            m_meshStructure,
                                                            numberOfDRCopyFaces,
                                                            numberOfDRInteriorFaces );

  delete[] numberOfDRCopyFaces;
  delete[] numberOfDRInteriorFaces;

  m_ltsTree = seissol::SeisSol::main.getMemoryManager().getLtsTree();
  m_lts = seissol::SeisSol::main.getMemoryManager().getLts();

  unsigned* ltsToMesh;
  unsigned numberOfMeshCells;
  // get cell information & mappings
  seissol::SeisSol::main.getLtsLayout().getCellInformation( m_ltsTree->var(m_lts->cellInformation),
                                                            ltsToMesh,
                                                            numberOfMeshCells );

  m_ltsLut.createLuts(  m_ltsTree,
                        ltsToMesh,
                        numberOfMeshCells );

  delete[] ltsToMesh;

  // derive lts setups
  seissol::initializers::time_stepping::deriveLtsSetups( m_timeStepping.numberOfLocalClusters,
                                                         m_meshStructure,
                                                         m_ltsTree->var(m_lts->cellInformation) );

  // initialize memory layout
  seissol::SeisSol::main.getMemoryManager().initializeMemoryLayout(enableFreeSurfaceIntegration);

  // add clusters
  seissol::SeisSol::main.timeManager().addClusters( m_timeStepping,
                                                    m_meshStructure,
                                                    seissol::SeisSol::main.getMemoryManager() );

  // get backward coupling
  m_globalData = seissol::SeisSol::main.getMemoryManager().getGlobalData();
}

#if defined(USE_NETCDF) && !defined(NETCDF_PASSIVE)
void seissol::Interoperability::setupNRFPointSources( char const* fileName )
{
  SeisSol::main.sourceTermManager().loadSourcesFromNRF(
    fileName,
    seissol::SeisSol::main.meshReader(),
    m_ltsTree,
    m_lts,
    &m_ltsLut,
    seissol::SeisSol::main.timeManager()
  );
}
#endif

void seissol::Interoperability::setupFSRMPointSources( double const* momentTensor,
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
                                                       double const* timeHistories )
{
  SeisSol::main.sourceTermManager().loadSourcesFromFSRM(
    momentTensor,
    velocityComponent,
    numberOfSources,
    centres,
    strikes,
    dips,
    rakes,
    onsets,
    areas,
    timestep,
    numberOfSamples,
    timeHistories,
    seissol::SeisSol::main.meshReader(),
    m_ltsTree,
    m_lts,
    &m_ltsLut,
    seissol::SeisSol::main.timeManager()
  );
}

void seissol::Interoperability::initializeModel(  char*   materialFileName,
                                                  int     anelasticity,
                                                  int     plasticity,
                                                  double* materialVal,
                                                  double* bulkFriction,
                                                  double* plastCo,
                                                  double* iniStress )
{
  auto nElements = seissol::SeisSol::main.meshReader().getElements().size();

  seissol::initializers::ParameterDB parameterDB;
  parameterDB.addParameter("rho",    materialVal);
  parameterDB.addParameter("mu",     materialVal + nElements);
  parameterDB.addParameter("lambda", materialVal + nElements*2);
  
  if (anelasticity == 1) {
    parameterDB.addParameter("Qp",  materialVal + nElements*3);
    parameterDB.addParameter("Qs",  materialVal + nElements*4);
  }

  if (plasticity == 1) {
    parameterDB.addParameter("bulkFriction", bulkFriction);
    parameterDB.addParameter("plastCo",      plastCo);
    parameterDB.addParameter("s_xx",         iniStress+0, 6);
    parameterDB.addParameter("s_yy",         iniStress+1, 6);
    parameterDB.addParameter("s_zz",         iniStress+2, 6);
    parameterDB.addParameter("s_xy",         iniStress+3, 6);
    parameterDB.addParameter("s_yz",         iniStress+4, 6);
    parameterDB.addParameter("s_xz",         iniStress+5, 6);
  }
  
  seissol::initializers::ElementBarycentreGenerator queryGen(seissol::SeisSol::main.meshReader());
  parameterDB.evaluateModel(std::string(materialFileName), queryGen);
}

void seissol::Interoperability::fitAttenuation( double rho,
                                                double mu,
                                                double lambda,
                                                double Qp,
                                                double Qs,
                                                seissol::model::Material& material )
{
  constexpr size_t numMaterialVals = 3 + 4*NUMBER_OF_RELAXATION_MECHANISMS;
  double materialFortran[numMaterialVals];
  f_interoperability_fitAttenuation(m_domain, rho, mu, lambda, Qp, Qs, materialFortran);
  seissol::model::setMaterial(materialFortran, numMaterialVals, &material);
}

void seissol::Interoperability::initializeFault( char*   modelFileName,
                                                 int     gpwise,
                                                 double* bndPoints,
                                                 int     numberOfBndPoints )
{
  seissol::initializers::ParameterDB parameterDB;
  for (auto const& kv : m_faultParameters) {
    parameterDB.addParameter(kv.first, kv.second);
  }
  
  if (gpwise != 0) {
    seissol::initializers::FaultGPGenerator queryGen( seissol::SeisSol::main.meshReader(),
                                                      reinterpret_cast<double(*)[2]>(bndPoints),
                                                      numberOfBndPoints );
    parameterDB.evaluateModel(std::string(modelFileName), queryGen);
  } else {
    seissol::initializers::FaultBarycentreGenerator queryGen( seissol::SeisSol::main.meshReader(),
                                                              numberOfBndPoints );
    parameterDB.evaluateModel(std::string(modelFileName), queryGen);
  }
}

void seissol::Interoperability::enableDynamicRupture() {
  // DR is always enabled if there are dynamic rupture cells
}

void seissol::Interoperability::setMaterial(int i_meshId, int i_side, double* i_materialVal, int i_numMaterialVals)
{
  int side = i_side - 1;
  seissol::model::Material* material;

  if (side < 0) {
    material = &m_ltsLut.lookup(m_lts->material, i_meshId - 1).local;
  } else {
    assert(side < 4);
    material = &m_ltsLut.lookup(m_lts->material, i_meshId - 1).neighbor[side];
  }

  seissol::model::setMaterial(i_materialVal, i_numMaterialVals, material);
}

#ifdef USE_PLASTICITY
void seissol::Interoperability::setInitialLoading( int i_meshId, double *i_initialLoading ) {
  PlasticityData& plasticity = m_ltsLut.lookup(m_lts->plasticity, i_meshId - 1);

  for( unsigned int l_stress = 0; l_stress < 6; l_stress++ ) {
    unsigned l_basis = 0;
    plasticity.initialLoading[l_stress] = i_initialLoading[ l_stress*NUMBER_OF_BASIS_FUNCTIONS + l_basis ];
  }
}
//synchronize element dependent plasticity parameters
void seissol::Interoperability::setPlasticParameters( int i_meshId, double* i_plasticParameters) {
  PlasticityData& plasticity = m_ltsLut.lookup(m_lts->plasticity, i_meshId - 1);
  CellMaterialData& material = m_ltsLut.lookup(m_lts->material, i_meshId - 1);

  double angularFriction = atan(i_plasticParameters[1]);

  plasticity.cohesionTimesCosAngularFriction = i_plasticParameters[0] * cos(angularFriction);
  plasticity.sinAngularFriction = sin(angularFriction);
  plasticity.mufactor = 1.0 / (2.0 * material.local.mu);

}


void seissol::Interoperability::setTv(double tv) {
  seissol::SeisSol::main.timeManager().setTv(tv);
}
#endif

void seissol::Interoperability::initializeCellLocalMatrices()
{
  // \todo Move this to some common initialization place
  seissol::initializers::initializeCellLocalMatrices( seissol::SeisSol::main.meshReader(),
                                                      m_ltsTree,
                                                      m_lts,
                                                      &m_ltsLut );

  seissol::initializers::initializeDynamicRuptureMatrices( seissol::SeisSol::main.meshReader(),
                                                           m_ltsTree,
                                                           m_lts,
                                                           &m_ltsLut,
                                                           seissol::SeisSol::main.getMemoryManager().getDynamicRuptureTree(),
                                                           seissol::SeisSol::main.getMemoryManager().getDynamicRupture(),
                                                           m_ltsFaceToMeshFace,
                                                           *seissol::SeisSol::main.getMemoryManager().getGlobalData(),
                                                           m_timeStepping );
}

template<typename T>
void seissol::Interoperability::synchronize(seissol::initializers::Variable<T> const& handle)
{
  unsigned *const (&meshToLts)[seissol::initializers::Lut::MaxDuplicates] = m_ltsLut.getMeshToLtsLut(handle.mask);
  unsigned* duplicatedMeshIds = m_ltsLut.getDuplicatedMeshIds(handle.mask);
  unsigned numberOfDuplicatedMeshIds = m_ltsLut.getNumberOfDuplicatedMeshIds(handle.mask);
  T* var = m_ltsTree->var(handle);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned dupMeshId = 0; dupMeshId < numberOfDuplicatedMeshIds; ++dupMeshId) {
    unsigned meshId = duplicatedMeshIds[dupMeshId];
    T* ref = &var[ meshToLts[0][meshId] ];
    for (unsigned dup = 1; dup < seissol::initializers::Lut::MaxDuplicates && meshToLts[dup][meshId] != std::numeric_limits<unsigned>::max(); ++dup) {
      memcpy(&var[ meshToLts[dup][meshId] ], ref, sizeof(T));
    }
  }
}

void seissol::Interoperability::synchronizeCellLocalData() {
  synchronize(m_lts->material);
#ifdef USE_PLASTICITY
  synchronize(m_lts->plasticity);
#endif
}

void seissol::Interoperability::synchronizeCopyLayerDofs() {
  synchronize(m_lts->dofs);
}

void seissol::Interoperability::enableWaveFieldOutput( double i_waveFieldInterval, const char *i_waveFieldFilename ) {
  seissol::SeisSol::main.waveFieldWriter().setWaveFieldInterval( i_waveFieldInterval );
  seissol::SeisSol::main.waveFieldWriter().enable();
  seissol::SeisSol::main.waveFieldWriter().setFilename( i_waveFieldFilename );
}

void seissol::Interoperability::enableFreeSurfaceOutput(int maxRefinementDepth)
{
	seissol::SeisSol::main.freeSurfaceWriter().enable();

	seissol::SeisSol::main.freeSurfaceIntegrator().initialize( maxRefinementDepth,
								m_lts,
								m_ltsTree,
								&m_ltsLut );
}


void seissol::Interoperability::enableCheckPointing( double i_checkPointInterval,
		const char *i_checkPointFilename, const char *i_checkPointBackend ) {
  seissol::SeisSol::main.simulator().setCheckPointInterval( i_checkPointInterval );
  if (strcmp(i_checkPointBackend, "posix") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::POSIX);
  else if (strcmp(i_checkPointBackend, "hdf5") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::HDF5);
  else if (strcmp(i_checkPointBackend, "mpio") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::MPIO);
  else if (strcmp(i_checkPointBackend, "mpio_async") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::MPIO_ASYNC);
  else if (strcmp(i_checkPointBackend, "sionlib") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::SIONLIB);
  else
	  logError() << "Unknown checkpoint backend";
  seissol::SeisSol::main.checkPointManager().setFilename( i_checkPointFilename );
}

void seissol::Interoperability::getIntegrationMask( int* i_integrationMask ) {
  seissol::SeisSol::main.postProcessor().setIntegrationMask(i_integrationMask);
}

void seissol::Interoperability::initializeIO(
		double* mu, double* slipRate1, double* slipRate2,
		double* slip, double* slip1, double* slip2, double* state, double* strength,
		int numSides, int numBndGP, int refinement, int* outputMask,
		double* outputRegionBounds,
		double freeSurfaceInterval, const char* freeSurfaceFilename,
    char const* xdmfWriterBackend,
    double receiverSamplingInterval, double receiverSyncInterval)
{
  auto type = writer::backendType(xdmfWriterBackend);
  
	// Initialize checkpointing
	int faultTimeStep;
	bool hasCheckpoint = seissol::SeisSol::main.checkPointManager().init(reinterpret_cast<real*>(m_ltsTree->var(m_lts->dofs)),
			m_ltsTree->getNumberOfCells(m_lts->dofs.mask) * tensor::Q::size(),
			mu, slipRate1, slipRate2, slip, slip1, slip2,
			state, strength, numSides, numBndGP,
			faultTimeStep);
	if (hasCheckpoint) {
		seissol::SeisSol::main.simulator().setCurrentTime(
			seissol::SeisSol::main.checkPointManager().header().time());
		seissol::SeisSol::main.faultWriter().setTimestep(faultTimeStep);
	}

  constexpr auto numberOfQuantities = tensor::Q::Shape[ sizeof(tensor::Q::Shape) / sizeof(tensor::Q::Shape[0]) - 1];

	// Initialize wave field output
	seissol::SeisSol::main.waveFieldWriter().init(
			numberOfQuantities, CONVERGENCE_ORDER,
			NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
			seissol::SeisSol::main.meshReader(),
			reinterpret_cast<const double*>(m_ltsTree->var(m_lts->dofs)),
			reinterpret_cast<const double*>(m_ltsTree->var(m_lts->pstrain)),
			seissol::SeisSol::main.postProcessor().getIntegrals(m_ltsTree),
			m_ltsLut.getMeshToLtsLut(m_lts->dofs.mask)[0],
			refinement, outputMask, outputRegionBounds,
      type);

	// Initialize free surface output
	seissol::SeisSol::main.freeSurfaceWriter().init(
		seissol::SeisSol::main.meshReader(),
		&seissol::SeisSol::main.freeSurfaceIntegrator(),
		freeSurfaceFilename, freeSurfaceInterval, type);

  auto& receiverWriter = seissol::SeisSol::main.receiverWriter();
  // Initialize receiver output
  receiverWriter.init(
    std::string(freeSurfaceFilename),
    receiverSamplingInterval,
    receiverSyncInterval
  );
  receiverWriter.addPoints(
    m_recPoints,
    seissol::SeisSol::main.meshReader(),
    m_ltsLut,
    *m_lts,
    m_globalData
  );
  seissol::SeisSol::main.timeManager().setReceiverClusters(receiverWriter);

	// I/O initialization is the last step that requires the mesh reader
	// (at least at the moment ...)

	// TODO(Lukas) Free the mesh reader if not doing convergence test.
	seissol::SeisSol::main.analysisWriter().init(&seissol::SeisSol::main.meshReader());
	//seissol::SeisSol::main.freeMeshReader();
}

void seissol::Interoperability::copyDynamicRuptureState()
{
	f_interoperability_copyDynamicRuptureState(m_domain);
}

void seissol::Interoperability::initInitialConditions()
{
  if (m_initialConditionType == "Planarwave") {
#ifdef MULTIPLE_SIMULATIONS
    for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
      m_iniConds.push_back(new physics::Planarwave((2.0*M_PI*s) / MULTIPLE_SIMULATIONS));
    }
#else
    m_iniConds.push_back(new physics::Planarwave());
#endif
  } else if (m_initialConditionType == "Zero") {
    m_iniConds.push_back(new physics::ZeroField());
  } else {
    throw std::runtime_error("Unknown initial condition type" + getInitialConditionType());
  }
}

void seissol::Interoperability::projectInitialField()
{
  initInitialConditions();

  if (m_initialConditionType == "Zero") {
    // Projection not necessary
    return;
  }
  initializers::projectInitialField(  getInitialConditions(),
                                      *m_globalData,
                                      seissol::SeisSol::main.meshReader(),
                                      *m_lts,
                                      m_ltsLut);
}

void seissol::Interoperability::getDofs( int    i_meshId,
                                         double o_dofs[tensor::QFortran::size()] ) {
  /// @yateto_todo: multiple sims?
  seissol::kernels::convertAlignedDofs( m_ltsLut.lookup(m_lts->dofs, i_meshId-1), o_dofs );
}

void seissol::Interoperability::getDofsFromDerivatives( int    i_meshId,
                                                        double o_dofs[tensor::QFortran::size()] ) {
  // assert that the cell provides derivatives
  assert( (m_ltsLut.lookup(m_lts->cellInformation, i_meshId-1).ltsSetup >> 9)%2 == 1 );

  // get DOFs from 0th derivatives
  seissol::kernels::convertAlignedDofs( m_ltsLut.lookup(m_lts->derivatives, i_meshId-1), o_dofs );
}

void seissol::Interoperability::getNeighborDofsFromDerivatives( int    i_meshId,
                                                                int    i_localFaceId,
                                                                double  o_dofs[tensor::QFortran::size()] ) {

  // get DOFs from 0th neighbors derivatives
  seissol::kernels::convertAlignedDofs(  m_ltsLut.lookup(m_lts->faceNeighbors, i_meshId-1)[ i_localFaceId-1 ],
                                         o_dofs );
}

seissol::initializers::Lut* seissol::Interoperability::getLtsLut() {
  return &m_ltsLut;
}

std::string seissol::Interoperability::getInitialConditionType() {
  return m_initialConditionType;
}

void seissol::Interoperability::simulate( double i_finalTime ) {
  seissol::SeisSol::main.simulator().setFinalTime( i_finalTime );

 seissol::SeisSol::main.simulator().simulate();
}

void seissol::Interoperability::finalizeIO()
{
	seissol::SeisSol::main.waveFieldWriter().close();
	seissol::SeisSol::main.checkPointManager().close();
	seissol::SeisSol::main.faultWriter().close();
	seissol::SeisSol::main.freeSurfaceWriter().close();
}

void seissol::Interoperability::faultOutput( double i_fullUpdateTime,
                                             double i_timeStepWidth )
{
  f_interoperability_faultOutput( m_domain, &i_fullUpdateTime, &i_timeStepWidth );
}

void seissol::Interoperability::evaluateFrictionLaw(  int face,
                                                      real godunov[CONVERGENCE_ORDER][seissol::tensor::godunovState::size()],
                                                      real imposedStatePlus[seissol::tensor::godunovState::size()],
                                                      real imposedStateMinus[seissol::tensor::godunovState::size()],
                                                      double i_fullUpdateTime,
                                                      double timePoints[CONVERGENCE_ORDER],
                                                      double timeWeights[CONVERGENCE_ORDER],
                                                      seissol::model::IsotropicWaveSpeeds const& waveSpeedsPlus,
                                                      seissol::model::IsotropicWaveSpeeds const& waveSpeedsMinus )
{
  int fFace = face + 1;
  int numberOfPoints = tensor::godunovState::Shape[0];
  int godunovLd = init::godunovState::Stop[0] - init::godunovState::Start[0];

  f_interoperability_evaluateFrictionLaw( m_domain,
                                          fFace,
                                         &godunov[0][0],
                                         &imposedStatePlus[0],
                                         &imposedStateMinus[0],
                                          numberOfPoints,
                                          godunovLd,
                                          &i_fullUpdateTime,
                                          &timePoints[0],
                                          &timeWeights[0],
                                          waveSpeedsPlus.density,
                                          waveSpeedsPlus.pWaveVelocity,
                                          waveSpeedsPlus.sWaveVelocity,
                                          waveSpeedsMinus.density,
                                          waveSpeedsMinus.pWaveVelocity,
                                          waveSpeedsMinus.sWaveVelocity);
}

void seissol::Interoperability::calcElementwiseFaultoutput(double time)
{
	f_interoperability_calcElementwiseFaultoutput(m_domain, time);
}


#ifdef USE_PLASTICITY
void seissol::Interoperability::computePlasticity(  double i_timeStep,
		                                            double *i_plasticParameters,
                                                    double (*i_initialLoading)[NUMBER_OF_BASIS_FUNCTIONS],
                                                    double *io_dofs,
													double *io_Energy,
													double *io_pstrain ) {
  // call fortran routine
  f_interoperability_computePlasticity(  m_domain,
                                        &i_timeStep,
										 NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
										 i_plasticParameters,
                                         i_initialLoading,
                                         io_dofs,
										 io_Energy,
										 io_pstrain );
}
#endif

void seissol::Interoperability::computeMInvJInvPhisAtSources(double x, double y, double z, unsigned element, real mInvJInvPhisAtSources[tensor::mInvJInvPhisAtSources::size()])
{
  double f_mInvJInvPhisAtSources[NUMBER_OF_BASIS_FUNCTIONS];

  int elem = static_cast<int>(element);
  f_interoperability_computeMInvJInvPhisAtSources(m_domain, x, y, z, elem, f_mInvJInvPhisAtSources);

  memset(mInvJInvPhisAtSources, 0, tensor::mInvJInvPhisAtSources::size() * sizeof(real));
  for (unsigned bf = 0; bf < NUMBER_OF_BASIS_FUNCTIONS; ++bf) {
    mInvJInvPhisAtSources[bf] = f_mInvJInvPhisAtSources[bf];
  }
}
