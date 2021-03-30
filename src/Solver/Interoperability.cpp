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

#include <cstddef>
#include <cstring>

#include "Interoperability.h"
#include "time_stepping/TimeManager.h"
#include "SeisSol.h"
#include <Initializer/CellLocalMatrices.h>
#include <Initializer/InitialFieldProjection.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/time_stepping/common.hpp>
#include <Initializer/typedefs.hpp>
#include <Equations/Setup.h>
#include <Numerical_aux/BasisFunction.h>
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

  void c_interoperability_initializeMemoryLayout(int clustering, bool enableFreeSurfaceIntegration) {
    e_interoperability.initializeMemoryLayout(clustering, enableFreeSurfaceIntegration);
  }

  void c_interoperability_initializeEasiBoundaries(char* fileName) {
    seissol::SeisSol::main.getMemoryManager().initializeEasiBoundaryReader(fileName);
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
                                            int     anisotropy,
                                            double* materialVal,
                                            double* bulkFriction,
                                            double* plastCo,
                                            double* iniStress,
                                            double* waveSpeeds )
  {
    e_interoperability.initializeModel( materialFileName,
                                        anelasticity,
                                        plasticity,
                                        anisotropy,
                                        materialVal,
                                        bulkFriction,
                                        plastCo,
                                        iniStress,
                                        waveSpeeds );
  }
  
  void c_interoperability_addFaultParameter(  char* name,
                                              double* memory  ) {
    e_interoperability.addFaultParameter(name, memory);
  }
  
  bool c_interoperability_faultParameterizedByTraction( char* modelFileName ) {
    return seissol::initializers::FaultParameterDB::faultParameterizedByTraction( std::string(modelFileName) );
  }

  bool c_interoperability_nucleationParameterizedByTraction( char* modelFileName ) {
    return seissol::initializers::FaultParameterDB::nucleationParameterizedByTraction( std::string(modelFileName) );
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

void c_interoperability_report_device_memory_status() {
  e_interoperability.reportDeviceMemoryStatus();
}

  void c_interoperability_deallocateMemoryManager() {
    e_interoperability.deallocateMemoryManager();
  }

  void c_interoperability_TetraDubinerP(double* phis, double xi, double eta, double zeta, int N) {
    assert(N > 0);
    unsigned idx = 0;
    for (unsigned int d = 0; d <= static_cast<unsigned>(N); ++d) {
      for (unsigned int k = 0; k <= d; ++k) {
        for (unsigned int j = 0; j <= d - k; ++j) {
            phis[idx++] = seissol::functions::TetraDubinerP({d - j - k, j, k}, {xi, eta, zeta});
        }
      }
    }
  }

  void c_interoperability_TriDubinerP(double* phis, double xi, double eta, int N) {
    assert(N > 0);
    unsigned idx = 0;
    for (unsigned int d = 0; d <= static_cast<unsigned>(N); ++d) {
      for (unsigned int j = 0; j <= d; ++j) {
        phis[idx++] = seissol::functions::TriDubinerP({d - j, j}, {xi, eta});
      }
    }
  }

  void c_interoperability_gradTriDubinerP(double* phis, double xi, double eta, int N) {
    assert(N > 0);
    unsigned idx = 0;
    for (unsigned int d = 0; d <= static_cast<unsigned>(N); ++d) {
      for (unsigned int j = 0; j <= d; ++j) {
        auto const grad = seissol::functions::gradTriDubinerP({d - j, j}, {xi, eta});
        for (auto const& g : grad) {
            phis[idx++] = g;
        }
      }
    }
  }

  double c_interoperability_M2invDiagonal(int no) {
      assert(no >= 0 && no < static_cast<int>(seissol::tensor::M2inv::Shape[0]));
      auto M2inv = seissol::init::M2inv::view::create(
        const_cast<real*>(seissol::init::M2inv::Values));
      return M2inv(no, no);
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
                                                      real*   i_QInterpolatedPlus,
                                                      real*   i_QInterpolatedMinus,
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
                                                      double  sWaveVelocityMinus,
                                                      real const* resampleMatrix );

  extern void f_interoperability_calcElementwiseFaultoutput( void *domain,
	                                                     double time );

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

  seissol::SeisSol::main.getMemoryManager().fixateLtsTree(m_timeStepping,
                                                          m_meshStructure,
                                                          numberOfDRCopyFaces,
                                                          numberOfDRInteriorFaces);

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


}

void seissol::Interoperability::initializeMemoryLayout(int clustering, bool enableFreeSurfaceIntegration) {
  // initialize memory layout
  seissol::SeisSol::main.getMemoryManager().initializeMemoryLayout(enableFreeSurfaceIntegration);

  // add clusters
  seissol::SeisSol::main.timeManager().addClusters( m_timeStepping,
                                                    m_meshStructure,
                                                    seissol::SeisSol::main.getMemoryManager() );

  // get backward coupling
  m_globalData = seissol::SeisSol::main.getMemoryManager().getGlobalDataOnHost();


  // initialize face lts trees
  seissol::SeisSol::main.getMemoryManager().fixateBoundaryLtsTree();
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
                                                  bool    anelasticity,
                                                  bool    plasticity,
                                                  bool    anisotropy,
                                                  double* materialVal,
                                                  double* bulkFriction,
                                                  double* plastCo,
                                                  double* iniStress,
                                                  double* waveSpeeds)
{
  //There are only some valid combinations of material properties
  // elastic materials
  // viscoelastic materials
  // elastoplastic materials
  // viscoplastic materials
  // anisotropic elastic materials
  

  //first initialize the (visco-)elastic part
  auto nElements = seissol::SeisSol::main.meshReader().getElements().size();
  seissol::initializers::ElementBarycentreGenerator queryGen(seissol::SeisSol::main.meshReader());
  auto calcWaveSpeeds = [&] (seissol::model::Material* material, int pos) {
    waveSpeeds[pos] = material->getMaxWaveSpeed();
    waveSpeeds[nElements + pos] = material->getSWaveSpeed();
    waveSpeeds[2*nElements + pos] = material->getSWaveSpeed();
  };
  if (anisotropy) { 
    if(anelasticity || plasticity) {
      logError() << "Anisotropy can not be combined with anelasticity or plasticity";
    }
    auto materials = std::vector<seissol::model::AnisotropicMaterial>(nElements);
    seissol::initializers::MaterialParameterDB<seissol::model::AnisotropicMaterial> parameterDB;
    parameterDB.setMaterialVector(&materials);
    parameterDB.evaluateModel(std::string(materialFileName), queryGen);
    for (unsigned int i = 0; i < nElements; i++) {
      materialVal[i] = materials[i].rho;
      materialVal[nElements + i] = materials[i].c11;
      materialVal[2*nElements + i] = materials[i].c12;
      materialVal[3*nElements + i] = materials[i].c13;
      materialVal[4*nElements + i] = materials[i].c14;
      materialVal[5*nElements + i] = materials[i].c15;
      materialVal[6*nElements + i] = materials[i].c16;
      materialVal[7*nElements + i] = materials[i].c22;
      materialVal[8*nElements + i] = materials[i].c23;
      materialVal[9*nElements + i] = materials[i].c24;
      materialVal[10*nElements + i] = materials[i].c25;
      materialVal[11*nElements + i] = materials[i].c26;
      materialVal[12*nElements + i] = materials[i].c33;
      materialVal[13*nElements + i] = materials[i].c34;
      materialVal[14*nElements + i] = materials[i].c35;
      materialVal[15*nElements + i] = materials[i].c36;
      materialVal[16*nElements + i] = materials[i].c44;
      materialVal[17*nElements + i] = materials[i].c45;
      materialVal[18*nElements + i] = materials[i].c46;
      materialVal[19*nElements + i] = materials[i].c55;
      materialVal[20*nElements + i] = materials[i].c56;
      materialVal[21*nElements + i] = materials[i].c66;
      calcWaveSpeeds(&materials[i], i);
    }
  } else {
    if (anelasticity) {
      auto materials = std::vector<seissol::model::ViscoElasticMaterial>(nElements);
      seissol::initializers::MaterialParameterDB<seissol::model::ViscoElasticMaterial> parameterDB;
      parameterDB.setMaterialVector(&materials);
      parameterDB.evaluateModel(std::string(materialFileName), queryGen);
      for (unsigned int i = 0; i < nElements; i++) {
        materialVal[i] = materials[i].rho;
        materialVal[nElements + i] = materials[i].mu;
        materialVal[2*nElements + i] = materials[i].lambda;
        materialVal[3*nElements + i] = materials[i].Qp;
        materialVal[4*nElements + i] = materials[i].Qs;
        calcWaveSpeeds(&materials[i], i);
      }
    } else {
      auto materials = std::vector<seissol::model::ElasticMaterial>(nElements);
      seissol::initializers::MaterialParameterDB<seissol::model::ElasticMaterial> parameterDB;
      parameterDB.setMaterialVector(&materials);
      parameterDB.evaluateModel(std::string(materialFileName), queryGen);
      for (unsigned int i = 0; i < nElements; i++) {
        materialVal[i] = materials[i].rho;
        materialVal[nElements + i] = materials[i].mu;
        materialVal[2*nElements + i] = materials[i].lambda;
        calcWaveSpeeds(&materials[i], i);
      }
    } 

    //now initialize the plasticity data
    if (plasticity) {
      auto materials = std::vector<seissol::model::Plasticity>(nElements);
      seissol::initializers::MaterialParameterDB<seissol::model::Plasticity> parameterDB;
      parameterDB.setMaterialVector(&materials);
      parameterDB.evaluateModel(std::string(materialFileName), queryGen);
      for (unsigned int i = 0; i < nElements; i++) {
        bulkFriction[i] = materials[i].bulkFriction;
        plastCo[i] = materials[i].plastCo;
        iniStress[i*6+0] = materials[i].s_xx;
        iniStress[i*6+1] = materials[i].s_yy;
        iniStress[i*6+2] = materials[i].s_zz;
        iniStress[i*6+3] = materials[i].s_xy;
        iniStress[i*6+4] = materials[i].s_yz;
        iniStress[i*6+5] = materials[i].s_xz;
      }
    } 
  }
}

void seissol::Interoperability::fitAttenuation( double rho,
                                                double mu,
                                                double lambda,
                                                double Qp,
                                                double Qs,
                                                seissol::model::ViscoElasticMaterial& material )
{
#if defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  constexpr size_t numMaterialVals = 3 + 4*NUMBER_OF_RELAXATION_MECHANISMS;
  double materialFortran[numMaterialVals];
  f_interoperability_fitAttenuation(m_domain, rho, mu, lambda, Qp, Qs, materialFortran);
  material =  seissol::model::ViscoElasticMaterial(materialFortran, numMaterialVals);
#endif
}

void seissol::Interoperability::initializeFault( char*   modelFileName,
                                                 int     gpwise,
                                                 double* bndPoints,
                                                 int     numberOfBndPoints )
{
  seissol::initializers::FaultParameterDB parameterDB;
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
  //TODO(Lukas)
  //Prepare for different materials in the same simulation
  //Use placement new because pointer to virtual function table gets overwritten by 0 during init.
#if defined USE_ANISOTROPIC
  new(material) seissol::model::AnisotropicMaterial(i_materialVal, i_numMaterialVals);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  new(material) seissol::model::ViscoElasticMaterial(i_materialVal, i_numMaterialVals);
#else 
  new(material) seissol::model::ElasticMaterial(i_materialVal, i_numMaterialVals);
#endif
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
#ifndef USE_ANISOTROPIC
  plasticity.mufactor = 1.0 / (2.0 * material.local.mu);
#else
  plasticity.mufactor = 3.0 / (2.0 * (material.local.c44 + material.local.c55 + material.local.c66));
#endif
}


void seissol::Interoperability::setTv(double tv) {
  seissol::SeisSol::main.timeManager().setTv(tv);
}
#endif

void seissol::Interoperability::initializeCellLocalMatrices()
{
  // \todo Move this to some common initialization place
  MeshReader& meshReader = seissol::SeisSol::main.meshReader();
  seissol::initializers::initializeCellLocalMatrices( meshReader,
                                                      m_ltsTree,
                                                      m_lts,
                                                      &m_ltsLut );

  initializers::MemoryManager& memoryManager = seissol::SeisSol::main.getMemoryManager();
  seissol::initializers::initializeDynamicRuptureMatrices( meshReader,
                                                           m_ltsTree,
                                                           m_lts,
                                                           &m_ltsLut,
                                                           memoryManager.getDynamicRuptureTree(),
                                                           memoryManager.getDynamicRupture(),
                                                           m_ltsFaceToMeshFace,
                                                           *memoryManager.getGlobalDataOnHost(),
                                                           m_timeStepping );

  seissol::initializers::initializeBoundaryMappings(meshReader,
                                                    memoryManager.getEasiBoundaryReader(),
                                                    m_ltsTree,
                                                    m_lts,
                                                    &m_ltsLut);

#ifdef ACL_DEVICE
  initializers::copyCellMatricesToDevice(m_ltsTree,
                                         m_lts,
                                         memoryManager.getDynamicRuptureTree(),
                                         memoryManager.getDynamicRupture(),
                                         memoryManager.getBoundaryTree(),
                                         memoryManager.getBoundary());

  memoryManager.recordExecutionPaths();
#endif
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
      memcpy(reinterpret_cast<void*>(&var[ meshToLts[dup][meshId] ]), ref, sizeof(T));
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
  if (kernels::size<tensor::Qane>() > 0) {
    synchronize(m_lts->dofsAne);
  }
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

  // record the clustering info i.e., distribution of elements within an LTS tree
  const std::vector<Element>& MeshElements = seissol::SeisSol::main.meshReader().getElements();
  std::vector<unsigned> LtsClusteringData(MeshElements.size());
  auto& LtsLayout = seissol::SeisSol::main.getLtsLayout();
  for (const auto& Element: MeshElements) {
    LtsClusteringData[Element.localId] = LtsLayout.getGlobalClusterId(Element.localId);
  }
	// Initialize wave field output
	seissol::SeisSol::main.waveFieldWriter().init(
      numberOfQuantities, CONVERGENCE_ORDER,
      NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
      seissol::SeisSol::main.meshReader(),
      LtsClusteringData,
      reinterpret_cast<const real*>(m_ltsTree->var(m_lts->dofs)),
      reinterpret_cast<const real*>(m_ltsTree->var(m_lts->pstrain)),
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
      m_iniConds.emplace_back(new physics::Planarwave(m_ltsLut.lookup(m_lts->material, 0), (2.0*M_PI*s) / MULTIPLE_SIMULATIONS));
    }
#else
    m_iniConds.emplace_back(new physics::Planarwave(m_ltsLut.lookup(m_lts->material, 0)));
#endif
  } else if (m_initialConditionType == "SuperimposedPlanarwave") {
#ifdef MULTIPLE_SIMULATIONS
    for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
      m_iniConds.emplace_back(new physics::SuperimposedPlanarwave(m_ltsLut.lookup(m_lts->material, 0), (2.0*M_PI*s) / MULTIPLE_SIMULATIONS));
    }
#else
    m_iniConds.emplace_back(new physics::SuperimposedPlanarwave(m_ltsLut.lookup(m_lts->material, 0)));
#endif
  } else if (m_initialConditionType == "Zero") {
    m_iniConds.emplace_back(new physics::ZeroField());
#if NUMBER_OF_RELAXATION_MECHANISMS == 0
  } else if (m_initialConditionType == "Scholte") {
    m_iniConds.emplace_back(new physics::ScholteWave());
  } else if (m_initialConditionType == "Snell") {
    m_iniConds.emplace_back(new physics::SnellsLaw());
  } else if (m_initialConditionType == "Ocean") {
    m_iniConds.emplace_back(new physics::Ocean());
#endif // NUMBER_OF_RELAXATION_MECHANISMS == 0
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

void seissol::Interoperability::deallocateMemoryManager() {
  seissol::SeisSol::main.deleteMemoryManager();
}

void seissol::Interoperability::faultOutput( double i_fullUpdateTime,
                                             double i_timeStepWidth )
{
  f_interoperability_faultOutput( m_domain, &i_fullUpdateTime, &i_timeStepWidth );
}

void seissol::Interoperability::evaluateFrictionLaw(  int face,
                                                      real QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()],
                                                      real QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()],
                                                      real imposedStatePlus[seissol::tensor::QInterpolated::size()],
                                                      real imposedStateMinus[seissol::tensor::QInterpolated::size()],
                                                      double i_fullUpdateTime,
                                                      double timePoints[CONVERGENCE_ORDER],
                                                      double timeWeights[CONVERGENCE_ORDER],
                                                      seissol::model::IsotropicWaveSpeeds const& waveSpeedsPlus,
                                                      seissol::model::IsotropicWaveSpeeds const& waveSpeedsMinus )
{
  int fFace = face + 1;
  int numberOfPoints = tensor::QInterpolated::Shape[0];
  int godunovLd = init::QInterpolated::Stop[0] - init::QInterpolated::Start[0];

  static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0], "Different number of quadrature points?");

  f_interoperability_evaluateFrictionLaw( m_domain,
                                          fFace,
                                         &QInterpolatedPlus[0][0],
                                         &QInterpolatedMinus[0][0],
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
                                          waveSpeedsMinus.sWaveVelocity,
                                          init::resample::Values );
}

void seissol::Interoperability::calcElementwiseFaultoutput(double time)
{
	f_interoperability_calcElementwiseFaultoutput(m_domain, time);
}

void seissol::Interoperability::reportDeviceMemoryStatus() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  constexpr size_t GB = 1024 * 1024 * 1024;
  const auto rank = seissol::MPI::mpi.rank();
  if (device.api->getCurrentlyOccupiedMem() > device.api->getMaxAvailableMem()) {
    std::stringstream stream;

    stream << "Device(" << rank << ")  memory is overloaded.\n"
           << "Totally allocated device memory, GB: " << device.api->getCurrentlyOccupiedMem() / GB << '\n'
           << "Allocated unified memory, GB: " << device.api->getCurrentlyOccupiedUnifiedMem() / GB << '\n'
           << "Memory capacity of device, GB: " << device.api->getMaxAvailableMem() / GB;

    logError() << stream.str();
  }
  else {
    double fraction = device.api->getCurrentlyOccupiedMem() / static_cast<double>(device.api->getMaxAvailableMem());
    logInfo() << "occupied memory on device(" << rank << "): " << fraction * 100.0 << "%";
  }
#endif
}
