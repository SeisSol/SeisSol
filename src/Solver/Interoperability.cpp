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

#include <cstddef>
#include <cstring>

#include "Interoperability.h"
#include "time_stepping/TimeManager.h"
#include "SeisSol.h"
#include <Physics/PointSource.h>
#include <Initializer/CellLocalMatrices.h>
#include <Model/Setup.h>

seissol::Interoperability e_interoperability;

template<typename T>
class index_sort_by_value
{
private:
    T const* value;
public:
    index_sort_by_value(T const* value) : value(value) {}
    inline bool operator()(unsigned i, unsigned j) const {
        return value[i] < value[j];
    }
};

/*
 * C bindings
 */
extern "C" {
  // fortran to c
  void c_interoperability_setDomain( void *i_domain ) {
    e_interoperability.setDomain( i_domain );
  }

  void c_interoperability_setTimeStepWidth( int    *i_meshId,
                                            double *i_timeStepWidth ) {
    e_interoperability.setTimeStepWidth( i_meshId,
                                         i_timeStepWidth );
  }

  void c_interoperability_initializeClusteredLts( int *i_clustering ) {
    e_interoperability.initializeClusteredLts( i_clustering );
  }
                                           
  void c_interopability_allocatePointSources(int* i_meshIds, int* i_numberOfPointSources)
  {
    e_interoperability.allocatePointSources(i_meshIds, i_numberOfPointSources);
  }
  
  void c_interopability_setupPointSource(int* i_source,
                                         double* i_mInvJInvPhisAtSources,
                                         double* i_localMomentTensor,
                                         double* i_strike,
                                         double* i_dip,
                                         double* i_rake,
                                         double* i_samples,
                                         int* i_numberOfSamples,
                                         double* i_onsetTime,
                                         double* i_samplingInterval )
  {
    e_interoperability.setupPointSource( i_source,
                                         i_mInvJInvPhisAtSources,
                                         i_localMomentTensor,
                                         i_strike,
                                         i_dip,
                                         i_rake,
                                         i_samples,
                                         i_numberOfSamples,
                                         i_onsetTime,
                                         i_samplingInterval );
  }

  void c_interoperability_addReceiver( int *i_receiverId,
                                       int *i_meshId ) {
    e_interoperability.addReceiver( i_receiverId,
                                    i_meshId );
  }

  void c_interoperability_setReceiverSampling( double *i_receiverSampling ) {
    e_interoperability.setReceiverSampling( i_receiverSampling );
  }

  void c_interoperability_enableDynamicRupture() {
    e_interoperability.enableDynamicRupture();
  }
  
  void c_interoperability_setMaterial( int*    i_meshId,
                                       int*    i_side,
                                       double* i_materialVal,
                                       int*    i_numMaterialVals ) {
    e_interoperability.setMaterial(i_meshId, i_side, i_materialVal, i_numMaterialVals);
  }

#ifdef USE_PLASTICITY
 void c_interoperability_setInitialLoading( int    *i_meshId,
                                            double *i_initialLoading ) {
    e_interoperability.setInitialLoading( i_meshId, i_initialLoading );
  }
#endif
  
  void c_interoperability_initializeCellLocalMatrices() {
    e_interoperability.initializeCellLocalMatrices();
  }

  void c_interoperability_synchronizeMaterial() {
    e_interoperability.synchronizeMaterial();
  }

  void c_interoperability_synchronizeCopyLayerDofs() {
    e_interoperability.synchronizeCopyLayerDofs();
  }

  void c_interoperability_enableWaveFieldOutput( double *i_waveFieldInterval, const char* i_waveFieldFilename ) {
    e_interoperability.enableWaveFieldOutput( i_waveFieldInterval, i_waveFieldFilename );
  }

  void c_interoperability_enableCheckPointing( double *i_checkPointInterval,
		  const char* i_checkPointFilename, const char* i_checkPointBackend ) {
    e_interoperability.enableCheckPointing( i_checkPointInterval,
    		i_checkPointFilename, i_checkPointBackend );
  }

  void c_interoperability_initializeIO( double* mu, double* slipRate1, double* slipRate2,
		  double* slip1, double* slip2, double* state, double* strength,
		  int *numSides, int *numBndGP) {
	  e_interoperability.initializeIO(mu, slipRate1, slipRate2, slip1, slip2, state, strength,
			  *numSides, *numBndGP);
  }

  void c_interoperability_addToDofs( int    *i_meshId,
                                     double  i_update[NUMBER_OF_DOFS] ) {
    e_interoperability.addToDofs( i_meshId, i_update );
  }

  void c_interoperability_getTimeDerivatives( int    *i_meshId,
                                              double  o_timeDerivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] ) {
    e_interoperability.getTimeDerivatives( i_meshId,
                                           o_timeDerivatives );
  }

  void c_interoperability_getFaceDerInt( int    *i_meshId,
                                         int    *i_localFaceId,
                                         double *i_timeStepWidth,
                                         double  o_timeDerivativesCell[CONVERGENCE_ORDER][NUMBER_OF_DOFS],
                                         double  o_timeDerivativesNeighbor[CONVERGENCE_ORDER][NUMBER_OF_DOFS],
                                         double  o_timeIntegratedCell[NUMBER_OF_DOFS],
                                         double  o_timeIntegratedNeighbor[NUMBER_OF_DOFS] ) {
    e_interoperability.getFaceDerInt( i_meshId,
                                      i_localFaceId,
                                      i_timeStepWidth,
                                      o_timeDerivativesCell,
                                      o_timeDerivativesNeighbor,
                                      o_timeIntegratedCell,
                                      o_timeIntegratedNeighbor );
  }

  void c_interoperability_getDofs( int    *i_meshId,
                                   double  o_timeDerivatives[NUMBER_OF_DOFS] ) {
    e_interoperability.getDofs( i_meshId, o_timeDerivatives );
  }

  void c_interoperability_getDofsFromDerivatives( int    *i_meshId,
                                                  double  o_dofs[NUMBER_OF_DOFS] ) {
    e_interoperability.getDofsFromDerivatives( i_meshId, o_dofs );
  }

  void c_interoperability_getNeighborDofsFromDerivatives( int    *i_meshId,
                                                          int    *i_localFaceId,
                                                          double  o_dofs[NUMBER_OF_DOFS] ) {
    e_interoperability.getNeighborDofsFromDerivatives( i_meshId, i_localFaceId, o_dofs );
  }

  void c_interoperability_simulate( double *i_finalTime ) {
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

  extern void f_interoperability_setDynamicRuptureTimeStep( void *i_domain,
		   	   	   	   	   	   	   	   	   	   	   	   	    int  *i_timeStep );

  extern void f_interoperability_getDynamicRuptureTimeStep( void *i_domain,
		  	  	  	  	  	  	  	  	  	  	  	  	    int  *i_timeStep );

  extern void f_interoperability_computeDynamicRupture( void   *i_domain,
                                                        double *i_fullUpdateTime,
                                                        double *i_timeStepWidth );

  extern void f_interoperability_computePlasticity( void    *i_domain,
                                                    double  *i_timestep,
                                                    double (*i_initialLoading)[NUMBER_OF_BASIS_FUNCTIONS],
                                                    double  *i_stresses,
                                                    double  *o_plasticUpdate );


  extern void f_interoperability_writeReceivers( void   *i_domain,
                                                 double *i_fullUpdateTime,
                                                 double *i_timeStepWidth,
                                                 double *i_receiverTime,
                                                 int    *i_numberOfReceivers, 
                                                 int    *i_receiverIds );
}

/*
 * C++ functions
 */
seissol::Interoperability::Interoperability() :
  m_domain(NULL), // reset domain pointer
  m_pointSourceToCluster(NULL),
  m_pointSources(NULL),
  m_cellToPointSources(NULL),
  m_numberOfCellToPointSourcesMappings(NULL)
{
}

seissol::Interoperability::~Interoperability()
{
  if (m_pointSources != NULL) {
    for (unsigned cluster = 0; cluster < m_timeStepping.numberOfLocalClusters; ++cluster) {
      for (unsigned pwlf = 0; pwlf < m_pointSources[cluster].numberOfSources; ++pwlf) {
        delete[] m_pointSources[cluster].momentTimeFunctions[pwlf].slopes;
        delete[] m_pointSources[cluster].momentTimeFunctions[pwlf].intercepts;
      }

      delete[] m_pointSources[cluster].mInvJInvPhisAtSources;
      delete[] m_pointSources[cluster].momentTensors;
      delete[] m_pointSources[cluster].momentTimeFunctions;  
    }
  }
  delete[] m_pointSources;
  delete[] m_pointSourceToCluster;
  
  if (m_cellToPointSources != NULL) {
    for (unsigned cluster = 0; cluster < m_timeStepping.numberOfLocalClusters; ++cluster) {
      delete[] m_cellToPointSources[cluster];
    }
  }
  delete[] m_cellToPointSources;
  delete[] m_numberOfCellToPointSourcesMappings;
}

void seissol::Interoperability::setDomain( void* i_domain ) {
  assert( i_domain != NULL );

  m_domain = i_domain;
}

void seissol::Interoperability::setTimeStepWidth( int    *i_meshId,
                                                  double *i_timeStepWidth ) {
  seissol::SeisSol::main.getLtsLayout().setTimeStepWidth( (*i_meshId)-1, *i_timeStepWidth );
}

void seissol::Interoperability::initializeClusteredLts( int *i_clustering ) {
  // assert a valid clustering
  assert( *i_clustering > 0 );

  // either derive a GTS or LTS layout
  if( *i_clustering == 1 ) {
    seissol::SeisSol::main.getLtsLayout().deriveLayout( single, 1);
  }
  else {
    seissol::SeisSol::main.getLtsLayout().deriveLayout( multiRate, *i_clustering );
  }

  // get cell information & mappings
  seissol::SeisSol::main.getLtsLayout().getCellInformation( m_numberOfMeshCells,
                                                            m_numberOfLtsCells,
                                                            m_numberOfCopyInteriorCells,
                                                            m_cellInformation,
                                                            m_meshToLts,
                                                            m_meshToCopyInterior,
                                                            m_meshToClusters,
                                                            m_ltsToMesh,
                                                            m_copyInteriorToMesh );

  // get the mesh structure
  seissol::SeisSol::main.getLtsLayout().getMeshStructure( m_meshStructure );

  // get time stepping
  seissol::SeisSol::main.getLtsLayout().getCrossClusterTimeStepping( m_timeStepping );

  // add clusters
  seissol::SeisSol::main.timeManager().addClusters( m_timeStepping,
                                                    m_meshStructure,
                                                    m_cellInformation,
                                                    m_meshToClusters );

  // get backward coupling
  seissol::SeisSol::main.timeManager().getRawData( m_globalData,
                                                   m_cellData,
                                                   m_dofs,
                                                   m_buffers,
                                                   m_derivatives,
                                                   m_faceNeighbors );
}

void seissol::Interoperability::allocatePointSources( int* i_meshIds,
                                                      int* i_numberOfPointSources )
{
  m_pointSourceToCluster = new unsigned[*i_numberOfPointSources][2];
  m_pointSources = new PointSources[m_timeStepping.numberOfLocalClusters];
  m_cellToPointSources = new CellToPointSourcesMapping*[m_timeStepping.numberOfLocalClusters];
  m_numberOfCellToPointSourcesMappings = new unsigned[m_timeStepping.numberOfLocalClusters];
  
  for (unsigned cluster = 0; cluster < m_timeStepping.numberOfLocalClusters; ++cluster) {
    m_pointSources[cluster].numberOfSources = 0;
    m_cellToPointSources[cluster] = NULL;
  }
  
  unsigned* sortedPointSourceIndex = new unsigned[*i_numberOfPointSources];
  for (unsigned source = 0; source < *i_numberOfPointSources; ++source) {
    sortedPointSourceIndex[source] = source;
  }
  std::sort(sortedPointSourceIndex, sortedPointSourceIndex + *i_numberOfPointSources, index_sort_by_value<int>(i_meshIds));
  
  std::vector< std::vector<unsigned> > clusterToPointSources(m_timeStepping.numberOfLocalClusters);

  // distribute sources to clusters
  for (unsigned source = 0; source < *i_numberOfPointSources; ++source) {
    unsigned sortedSource = sortedPointSourceIndex[source];
    // get cell id (Fortran-formatting expected)
    unsigned meshId = i_meshIds[sortedSource]-1;
    unsigned cluster = m_meshToClusters[meshId][0];
    m_pointSourceToCluster[sortedSource][0] = cluster;
    m_pointSourceToCluster[sortedSource][1] = m_pointSources[cluster].numberOfSources++;
    clusterToPointSources[cluster].push_back(sortedSource);
  }
  
  delete[] sortedPointSourceIndex;

  unsigned clusterOffset = 0;
  for (unsigned cluster = 0; cluster < m_timeStepping.numberOfLocalClusters; ++cluster) {
    unsigned numberOfSourcesInCluster = m_pointSources[cluster].numberOfSources;
    
    // Find the cell offsets for a point source. As a cell has 4 neighbors,
    // the cell might exist up to 4 times in the copy layer.
    CellToPointSourcesMapping* cellToPointSources = new CellToPointSourcesMapping[4 * numberOfSourcesInCluster + 1];

    int mapping = -1;
    unsigned lastMeshId = std::numeric_limits<unsigned>::max();
    // add only the interior layer offsets
    for (std::vector<unsigned>::const_iterator source = clusterToPointSources[cluster].begin(); source != clusterToPointSources[cluster].end(); ++source) {
      unsigned meshId = i_meshIds[*source]-1;      
      unsigned offset = m_meshToCopyInterior[meshId];
      // If we have a interior cell
      if (offset >= clusterOffset + m_meshStructure[cluster].numberOfCopyCells) {
        unsigned clusterPointSourceId = m_pointSourceToCluster[*source][1];
        if (lastMeshId == meshId) {
          assert(clusterPointSourceId <= cellToPointSources[mapping].pointSourcesOffset + cellToPointSources[mapping].numberOfPointSources);
          ++cellToPointSources[mapping].numberOfPointSources;
        } else {
          lastMeshId = meshId;
          ++mapping;
          cellToPointSources[mapping].copyInteriorOffset = offset - clusterOffset;
          cellToPointSources[mapping].pointSourcesOffset = clusterPointSourceId;
          cellToPointSources[mapping].numberOfPointSources = 1;
        }
      }
    }
    
    // add the copy layer offsets
    for (unsigned cell = 0; cell < m_meshStructure[cluster].numberOfCopyCells; ++cell) {
      unsigned cellMeshId = m_copyInteriorToMesh[cell + clusterOffset];
      assert(mapping < 4 * static_cast<int>(numberOfSourcesInCluster));
      ++mapping;
      cellToPointSources[mapping].numberOfPointSources = 0;
      cellToPointSources[mapping].copyInteriorOffset = cell;
      
      for (std::vector<unsigned>::const_iterator source = clusterToPointSources[cluster].begin(); source != clusterToPointSources[cluster].end(); ++source) {
        // get cell id (Fortran-formatting expected)
        unsigned sourceMeshId = i_meshIds[*source]-1;
        if (sourceMeshId == cellMeshId) {
          unsigned clusterPointSourceId = m_pointSourceToCluster[*source][1];
          if (cellToPointSources[mapping].numberOfPointSources == 0) {
            cellToPointSources[mapping].pointSourcesOffset = clusterPointSourceId;
          }
          assert(clusterPointSourceId <= cellToPointSources[mapping].pointSourcesOffset + cellToPointSources[mapping].numberOfPointSources);
          ++cellToPointSources[mapping].numberOfPointSources;
        }
      }
        
      if (cellToPointSources[mapping].numberOfPointSources == 0) {
        --mapping;
      }
    }
    
    m_numberOfCellToPointSourcesMappings[cluster] = mapping+1;
    
    m_cellToPointSources[cluster] = new CellToPointSourcesMapping[ m_numberOfCellToPointSourcesMappings[cluster] ];
    for (unsigned i = 0; i < m_numberOfCellToPointSourcesMappings[cluster]; ++i) {
      m_cellToPointSources[cluster][i] = cellToPointSources[i];
    }    
    delete[] cellToPointSources;
    clusterOffset += m_meshStructure[cluster].numberOfCopyCells + m_meshStructure[cluster].numberOfInteriorCells;

    m_pointSources[cluster].mInvJInvPhisAtSources = new real[numberOfSourcesInCluster][NUMBER_OF_ALIGNED_BASIS_FUNCTIONS];
    m_pointSources[cluster].momentTensors = new real[numberOfSourcesInCluster][NUMBER_OF_QUANTITIES];
    m_pointSources[cluster].momentTimeFunctions = new PiecewiseLinearFunction1D[numberOfSourcesInCluster];
    
    for (unsigned pwlf = 0; pwlf < numberOfSourcesInCluster; ++pwlf) {
      m_pointSources[cluster].momentTimeFunctions[pwlf].slopes = NULL;
      m_pointSources[cluster].momentTimeFunctions[pwlf].intercepts = NULL;
    }
  }
  
  seissol::SeisSol::main.timeManager().setPointSourcesForClusters( m_cellToPointSources,
                                                                   m_numberOfCellToPointSourcesMappings,
                                                                   m_pointSources,
                                                                   m_timeStepping.numberOfLocalClusters );
}

void seissol::Interoperability::setupPointSource( int* i_source,
                                                  double* i_mInvJInvPhisAtSources,
                                                  double* i_localMomentTensor,
                                                  double* i_strike,
                                                  double* i_dip,
                                                  double* i_rake,
                                                  double* i_samples,
                                                  int* i_numberOfSamples,
                                                  double* i_onsetTime,
                                                  double* i_samplingInterval )
{  
  // get point source id (Fortan-formatting expected)
  unsigned cluster = m_pointSourceToCluster[*i_source-1][0];
  unsigned localSourceId = m_pointSourceToCluster[*i_source-1][1];
  
  memset(m_pointSources[cluster].mInvJInvPhisAtSources[localSourceId], 0.0, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * sizeof(real));
  seissol::kernels::copySubMatrix( i_mInvJInvPhisAtSources, NUMBER_OF_BASIS_FUNCTIONS, 1, NUMBER_OF_BASIS_FUNCTIONS,
                                   m_pointSources[cluster].mInvJInvPhisAtSources[localSourceId], NUMBER_OF_BASIS_FUNCTIONS, 1, NUMBER_OF_BASIS_FUNCTIONS );

  // Do floating point precision conversion if necessary
  real localMomentTensor[3][3]; 
  seissol::kernels::copySubMatrix( i_localMomentTensor, 3, 3, 3,
                                     &localMomentTensor[0][0], 3, 3, 3 );

  seissol::physics::transformMomentTensor( localMomentTensor,
                                           *i_strike,
                                           *i_dip,
                                           *i_rake,
                                           m_pointSources[cluster].momentTensors[localSourceId] );  

  seissol::physics::samplesToPiecewiseLinearFunction1D( i_samples,
                                                        *i_numberOfSamples,
                                                        *i_onsetTime,
                                                        *i_samplingInterval,
                                                        &m_pointSources[cluster].momentTimeFunctions[localSourceId] );
}

void seissol::Interoperability::addReceiver( int *i_receiverId,
                                             int *i_meshId ) {
  assert( i_meshId     != NULL );
  assert( i_receiverId != NULL );

  seissol::SeisSol::main.timeManager().addReceiver( *i_receiverId, *i_meshId );
}

void seissol::Interoperability::setReceiverSampling( double *i_receiverSampling ) {
  assert( i_receiverSampling != NULL );
  assert( *i_receiverSampling > 0 );

  seissol::SeisSol::main.timeManager().setReceiverSampling( *i_receiverSampling );
}

void seissol::Interoperability::enableDynamicRupture() {
  seissol::SeisSol::main.timeManager().enableDynamicRupture();
}

void seissol::Interoperability::setMaterial(int* i_meshId, int* i_side, double* i_materialVal, int* i_numMaterialVals)
{
  unsigned int copyInteriorId = m_meshToCopyInterior[*i_meshId - 1];
  int side = *i_side - 1;
  seissol::model::Material* material;
  
  if (side < 0) {
    material = &m_cellData->material[copyInteriorId].local;
  } else {
    assert(side < 4);
    material = &m_cellData->material[copyInteriorId].neighbor[side];
  }
  
  seissol::model::setMaterial(i_materialVal, *i_numMaterialVals, material);
}

#ifdef USE_PLASTICITY
void seissol::Interoperability::setInitialLoading( int* i_meshId, double *i_initialLoading ) {\
  unsigned int l_copyInteriorId = m_meshToCopyInterior[*i_meshId - 1];

  for( unsigned int l_stress = 0; l_stress < 6; l_stress++ ) {
    for( unsigned int l_basis = 0; l_basis < NUMBER_OF_BASIS_FUNCTIONS; l_basis++ ) {
      m_cellData->neighboringIntegration[l_copyInteriorId].initialLoading[l_stress][l_basis] = i_initialLoading[ l_stress*NUMBER_OF_BASIS_FUNCTIONS + l_basis ];
    }
  }
}
#endif

void seissol::Interoperability::initializeCellLocalMatrices()
{
  // \todo Move this to some common initialization place
  seissol::initializers::initializeCellLocalMatrices( seissol::SeisSol::main.meshReader(),
                                                      m_copyInteriorToMesh,
                                                      m_meshToLts,
                                                      m_numberOfCopyInteriorCells,
                                                      m_cellInformation,
                                                      m_cellData );
}

void seissol::Interoperability::synchronizeMaterial() {
  // iterate over the mesh and set all redundant data
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < m_numberOfCopyInteriorCells; l_cell++ ) {
    unsigned sourceId = m_meshToCopyInterior[m_copyInteriorToMesh[l_cell]];
    m_cellData->material[l_cell].local = m_cellData->material[sourceId].local;
    for (unsigned side = 0; side < 4; ++side) {
      m_cellData->material[l_cell].neighbor[side] = m_cellData->material[sourceId].neighbor[side];
    }
  }
}

void seissol::Interoperability::synchronizeCopyLayerDofs() {
  unsigned int l_offset = 0;

  for( unsigned int l_cluster = 0; l_cluster < m_timeStepping.numberOfLocalClusters; l_cluster++ ) {
    // sync DOFs
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for( unsigned int l_cell = 0; l_cell < m_meshStructure[l_cluster].numberOfCopyCells; l_cell++ ) {
      for( int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
        m_dofs[l_cell+l_offset][l_dof] = m_dofs[ m_meshToCopyInterior[m_copyInteriorToMesh[l_cell+l_offset]] ][l_dof];
      }
    }

    // update offset
    l_offset += m_meshStructure[l_cluster].numberOfCopyCells + m_meshStructure[l_cluster].numberOfInteriorCells;
  }
}

void seissol::Interoperability::enableWaveFieldOutput( double *i_waveFieldInterval, const char *i_waveFieldFilename ) {
  seissol::SeisSol::main.simulator().setWaveFieldInterval( *i_waveFieldInterval );
  seissol::SeisSol::main.waveFieldWriter().enable();
  seissol::SeisSol::main.waveFieldWriter().setFilename( i_waveFieldFilename );
}

void seissol::Interoperability::enableCheckPointing( double *i_checkPointInterval,
		const char *i_checkPointFilename, const char *i_checkPointBackend ) {
  seissol::SeisSol::main.simulator().setCheckPointInterval( *i_checkPointInterval );
  if (strcmp(i_checkPointBackend, "posix") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::POSIX);
  else if (strcmp(i_checkPointBackend, "hdf5") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::HDF5);
  else if (strcmp(i_checkPointBackend, "mpio") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::MPIO);
  else if (strcmp(i_checkPointBackend, "mpio_async") == 0)
	  seissol::SeisSol::main.checkPointManager().setBackend(checkpoint::MPIO_ASYNC);
  else
	  logError() << "Unknown checkpoint backend";
  seissol::SeisSol::main.checkPointManager().setFilename( i_checkPointFilename );
}

void seissol::Interoperability::initializeIO(
		  double* mu, double* slipRate1, double* slipRate2,
		  double* slip1, double* slip2, double* state, double* strength,
		  int numSides, int numBndGP)
{
	  // Initialize checkpointing
	  double currentTime;
	  int waveFieldTimeStep = 0;
	  int faultTimeStep;
	  bool hasCheckpoint = seissol::SeisSol::main.checkPointManager().init(reinterpret_cast<real*>(m_dofs),
				  m_numberOfCopyInteriorCells * NUMBER_OF_ALIGNED_DOFS,
				  mu, slipRate1, slipRate2, slip1, slip2,
				  state, strength, numSides, numBndGP,
				  currentTime, waveFieldTimeStep, faultTimeStep);
	  if (hasCheckpoint) {
		  seissol::SeisSol::main.simulator().setCurrentTime(currentTime);
		  f_interoperability_setDynamicRuptureTimeStep(m_domain, &faultTimeStep);
	  }

	  // Initialize wave field output
	  seissol::SeisSol::main.waveFieldWriter().init(
			  NUMBER_OF_QUANTITIES, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
			  seissol::SeisSol::main.meshReader(),
			  reinterpret_cast<const double*>(m_dofs), m_meshToCopyInterior,
			  waveFieldTimeStep);

	  // I/O initialization is the last step that requires the mesh reader
	  // (at least at the moment ...)
	  seissol::SeisSol::main.freeMeshReader();
}

void seissol::Interoperability::getDynamicRuptureTimeStep(int &o_timeStep)
{
	f_interoperability_getDynamicRuptureTimeStep(m_domain, &o_timeStep);
}

void seissol::Interoperability::addToDofs( int    *i_meshId,
                                           double  i_update[NUMBER_OF_DOFS] ) {
  seissol::kernels::addToAlignedDofs( i_update, m_dofs[ m_meshToCopyInterior[(*i_meshId)-1] ] );
}

void seissol::Interoperability::getTimeDerivatives( int    *i_meshId,
                                                    double  o_timeDerivatives[CONVERGENCE_ORDER][NUMBER_OF_DOFS] ) {
  real l_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));
  real l_timeDerivatives[NUMBER_OF_ALIGNED_DERS] __attribute__((aligned(ALIGNMENT)));

  m_timeKernel.computeAder( 0,
                            m_globalData->stiffnessMatricesTransposed,
                            m_dofs[ m_meshToCopyInterior[(*i_meshId)-1] ],
                            m_cellData->localIntegration[ m_meshToCopyInterior[ (*i_meshId)-1] ].starMatrices,
#ifdef REQUIRE_SOURCE_MATRIX
                            m_cellData->localIntegration[ m_meshToCopyInterior[ (*i_meshId)-1] ].sourceMatrix,
#endif
                            l_timeIntegrated,
                            l_timeDerivatives );

// We cannot use derivative compression with a source matrix
#ifdef REQUIRE_SOURCE_MATRIX
  for (unsigned order = 0; order < CONVERGENCE_ORDER; ++order) {
    seissol::kernels::copySubMatrix( &l_timeDerivatives[order * NUMBER_OF_ALIGNED_DOFS],
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                                     o_timeDerivatives[order],
                                     NUMBER_OF_BASIS_FUNCTIONS,
                                     NUMBER_OF_QUANTITIES,
                                     NUMBER_OF_BASIS_FUNCTIONS );
  }
#else
  seissol::kernels::convertAlignedCompressedTimeDerivatives( l_timeDerivatives,
                                                             o_timeDerivatives );
#endif
}

void seissol::Interoperability::getFaceDerInt( int    *i_meshId,
                                               int    *i_localFaceId,
                                               double *i_timeStepWidth,
                                               double  o_timeDerivativesCell[CONVERGENCE_ORDER][NUMBER_OF_DOFS],
                                               double  o_timeDerivativesNeighbor[CONVERGENCE_ORDER][NUMBER_OF_DOFS],
                                               double  o_timeIntegratedCell[NUMBER_OF_DOFS],
                                               double  o_timeIntegratedNeighbor[NUMBER_OF_DOFS] ) {
  // assert that the cell provides derivatives
  assert( (m_cellInformation[ m_meshToLts[ (*i_meshId)-1 ] ].ltsSetup >> 9)%2 == 1 );

  // get cells derivatives
  seissol::kernels::convertAlignedCompressedTimeDerivatives( m_derivatives[ m_meshToLts[ (*i_meshId)-1 ] ],
                                                             o_timeDerivativesCell );

  // get neighbors derivatives
  seissol::kernels::convertAlignedCompressedTimeDerivatives( m_faceNeighbors[ m_meshToCopyInterior[ (*i_meshId)-1 ] ][ (*i_localFaceId)-1 ],
                                                             o_timeDerivativesNeighbor );

  real l_timeIntegrated[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));

  // compute time integrated DOFs of the cell
  m_timeKernel.computeIntegral(  0,
                                 0,
                                *i_timeStepWidth,
                                 m_derivatives[ m_meshToLts[ (*i_meshId)-1 ] ],
                                 l_timeIntegrated );

  seissol::kernels::convertAlignedDofs( l_timeIntegrated, o_timeIntegratedCell );

  // compute time integrated dofs of the neighbor

  m_timeKernel.computeIntegral(  0,
                                 0,
                                *i_timeStepWidth,
                                 m_faceNeighbors[ m_meshToCopyInterior[ (*i_meshId)-1 ] ][ (*i_localFaceId)-1 ],
                                 l_timeIntegrated );

  seissol::kernels::convertAlignedDofs( l_timeIntegrated, o_timeIntegratedNeighbor );
}

void seissol::Interoperability::getDofs( int    *i_meshId,
                                         double  o_dofs[NUMBER_OF_DOFS] ) {
  unsigned int l_cellId = m_meshToCopyInterior[ (*i_meshId)-1 ];

  seissol::kernels::convertAlignedDofs( m_dofs[l_cellId], o_dofs );
}

void seissol::Interoperability::getDofsFromDerivatives( int    *i_meshId,
                                                        double  o_dofs[NUMBER_OF_DOFS] ) {
  // assert that the cell provides derivatives
  assert( (m_cellInformation[ m_meshToLts[ (*i_meshId)-1 ] ].ltsSetup >> 9)%2 == 1 );

  // get DOFs from 0th derivatives
  seissol::kernels::convertAlignedDofs( m_derivatives[ m_meshToLts[ (*i_meshId)-1 ] ], o_dofs );
}

void seissol::Interoperability::getNeighborDofsFromDerivatives( int    *i_meshId,
                                                                int    *i_localFaceId,
                                                                double  o_dofs[NUMBER_OF_DOFS] ) {

  // get DOFs from 0th neighbors derivatives
  seissol::kernels::convertAlignedDofs(  m_faceNeighbors[ m_meshToCopyInterior[ (*i_meshId)-1 ] ][ (*i_localFaceId)-1 ],
                                         o_dofs );
}

void seissol::Interoperability::simulate( double *i_finalTime ) {
  seissol::SeisSol::main.simulator().setFinalTime( *i_finalTime );

 seissol::SeisSol::main.simulator().simulate();
}

void seissol::Interoperability::finalizeIO()
{
	seissol::SeisSol::main.waveFieldWriter().close();
	seissol::SeisSol::main.checkPointManager().close();
}

void seissol::Interoperability::writeReceivers( double i_fullUpdateTime,
                                                double i_timeStepWidth,
                                                double i_receiverTime,
                                                std::vector< int > &i_receiverIds ) {
  assert( i_receiverIds.size() > 0 );

  // get number of elements
  int l_numberOfReceivers = i_receiverIds.size();

  // get a pointer to the raw data
  int *l_receiverIds = &i_receiverIds.front();

  // call the fortran routine
  f_interoperability_writeReceivers( m_domain,
                                    &i_fullUpdateTime,
                                    &i_timeStepWidth,
                                    &i_receiverTime,
                                    &l_numberOfReceivers, 
                                     l_receiverIds );
}

void seissol::Interoperability::computeDynamicRupture( double i_fullUpdateTime,
                                                       double i_timeStepWidth ) {
  f_interoperability_computeDynamicRupture(  m_domain,
                                            &i_fullUpdateTime,
                                            &i_timeStepWidth );
}

#ifdef USE_PLASTICITY
void seissol::Interoperability::computePlasticity(  double i_timeStep,
                                                    double (*i_initialLoading)[NUMBER_OF_BASIS_FUNCTIONS],
                                                    double *io_dofs ) {
  // collect stresses
  double l_stresses[6*NUMBER_OF_BASIS_FUNCTIONS];

  for( unsigned int l_quantity = 0; l_quantity < 6; l_quantity++ ) {
    for( unsigned int l_dof = 0; l_dof < NUMBER_OF_BASIS_FUNCTIONS; l_dof++ ) {
      l_stresses[l_quantity*NUMBER_OF_BASIS_FUNCTIONS+l_dof] = io_dofs[l_quantity * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + l_dof];
    }
  }

  // platic update
  double l_plasticUpdate[6*NUMBER_OF_BASIS_FUNCTIONS];

  // call fortran routine
  f_interoperability_computePlasticity(  m_domain,
                                        &i_timeStep,
                                         i_initialLoading,
                                         l_stresses,
                                         l_plasticUpdate );

  // update degrees of freedom
  for( unsigned int l_quantity = 0; l_quantity < 6; l_quantity++ ) {
    for( unsigned int l_dof = 0; l_dof < NUMBER_OF_BASIS_FUNCTIONS; l_dof++ ) {
      io_dofs[l_quantity * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS + l_dof] += l_plasticUpdate[l_quantity*NUMBER_OF_BASIS_FUNCTIONS + l_dof];
    }
  }
}
#endif
