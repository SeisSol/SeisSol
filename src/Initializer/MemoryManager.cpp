/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Memory management of SeisSol.
 **/
#include "SeisSol.h"
#include "MemoryManager.h"
#include "InternalState.h"

#include <Kernels/common.hpp>

#include <generated_code/init.h>

#ifdef _OPENMP
#include <omp.h>
#endif

seissol::initializers::MemoryManager::MemoryManager()
  : m_integrationBufferLTS(NULL)
{
}

void seissol::initializers::MemoryManager::initialize()
{
  // allocate thread-local LTS integration buffers
  allocateIntegrationBufferLTS();

  // initialize global matrices
#ifndef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
  initializeGlobalData( m_globalData );
#else
  // determine the number of threads
  unsigned int l_numberOfThreads = omp_get_max_threads();
  unsigned int l_numberOfCopiesCeil = (l_numberOfThreads%NUMBER_OF_THREADS_PER_GLOBALDATA_COPY == 0) ? 0 : 1;
  unsigned int l_numberOfCopies = (l_numberOfThreads/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY) + l_numberOfCopiesCeil;
  logInfo(0) << "Number of GlobalData copies: " << l_numberOfCopies;

  m_globalDataCopies = new GlobalData[l_numberOfCopies]; 

  // initialize in parallel to obtain best possible NUMA placement
  #pragma omp parallel
  {
    if (omp_get_thread_num()%NUMBER_OF_THREADS_PER_GLOBALDATA_COPY == 0) {
      // @TODO check why initializeGlobalMatrices is not thread-safe
      #pragma omp critical
      {
        initializeGlobalData( m_globalDataCopies[omp_get_thread_num()/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY] );
      }
    }
  }

  // set master structure
  m_globalData = m_globalDataCopies[0];
#endif
}

seissol::initializers::MemoryManager::~MemoryManager() {
  // free members
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
  delete[] m_globalDataCopies;
#endif
}

void seissol::initializers::MemoryManager::initializeGlobalData( struct GlobalData &o_globalData )
{
  // DG global matrices
  real* globalMatrixMem = static_cast<real*>(m_memoryAllocator.allocateMemory( seissol::model::globalMatrixOffsets[seissol::model::numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, MEMKIND_GLOBAL ));
  for (unsigned matrix = 0; matrix < seissol::model::numGlobalMatrices; ++matrix) {
    memcpy(
      &globalMatrixMem[ seissol::model::globalMatrixOffsets[matrix] ],
      seissol::model::globalMatrixValues[matrix],
      (seissol::model::globalMatrixOffsets[matrix+1] - seissol::model::globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    o_globalData.stiffnessMatricesTransposed[transposedStiffness] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[transposedStiffness] ];
  }
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    o_globalData.stiffnessMatrices[stiffness] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[3 + stiffness] ];
  }
  for (unsigned cob = 0; cob < 4; ++cob) {
    o_globalData.changeOfBasisMatrices[cob] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[6 + cob] ];
  }
  for (unsigned cob = 0; cob < 4; ++cob) {
    o_globalData.neighbourChangeOfBasisMatricesTransposed[cob] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[10 + cob] ];
  }
  for (unsigned cob = 0; cob < 4; ++cob) {
    o_globalData.localChangeOfBasisMatricesTransposed[cob] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[14 + cob] ];
  }
  for (unsigned flux = 0; flux < 3; ++flux) {
    o_globalData.neighbourFluxMatrices[flux] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[18 + flux] ];
  }

  // @TODO Integrate this step into the code generator
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = &globalMatrixMem[ seissol::model::globalMatrixOffsets[transposedStiffness] ];
    for (unsigned i = 0; i < seissol::model::globalMatrixOffsets[transposedStiffness+1]-seissol::model::globalMatrixOffsets[transposedStiffness]; ++i) {
      matrix[i] *= -1.0;
    }
  }

  o_globalData.integrationBufferLTS = m_integrationBufferLTS;

  // Dynamic Rupture global matrices  
  real* drGlobalMatrixMem = static_cast<real*>(m_memoryAllocator.allocateMemory( seissol::model::dr_globalMatrixOffsets[seissol::model::dr_numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, MEMKIND_GLOBAL ));
  for (unsigned matrix = 0; matrix < seissol::model::dr_numGlobalMatrices; ++matrix) {
    memcpy(
      &drGlobalMatrixMem[ seissol::model::dr_globalMatrixOffsets[matrix] ],
      seissol::model::dr_globalMatrixValues[matrix],
      (seissol::model::dr_globalMatrixOffsets[matrix+1] - seissol::model::dr_globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  for (unsigned face = 0; face < 4; ++face) {
    for (unsigned h = 0; h < 4; ++h) {
      o_globalData.nodalFluxMatrices[face][h] = &drGlobalMatrixMem[ seissol::model::dr_globalMatrixOffsets[4*face+h] ];
      o_globalData.faceToNodalMatrices[face][h] = &drGlobalMatrixMem[ seissol::model::dr_globalMatrixOffsets[16 + 4*face+h] ];
    }
  }
  
  real* plasticityGlobalMatrixMem = static_cast<real*>(m_memoryAllocator.allocateMemory( seissol::model::plasticity_globalMatrixOffsets[seissol::model::plasticity_numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, MEMKIND_GLOBAL ));
  for (unsigned matrix = 0; matrix < seissol::model::plasticity_numGlobalMatrices; ++matrix) {
    memcpy(
      &plasticityGlobalMatrixMem[ seissol::model::plasticity_globalMatrixOffsets[matrix] ],
      seissol::model::plasticity_globalMatrixValues[matrix],
      (seissol::model::plasticity_globalMatrixOffsets[matrix+1] - seissol::model::plasticity_globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  o_globalData.vandermondeMatrix = &plasticityGlobalMatrixMem[ seissol::model::plasticity_globalMatrixOffsets[0] ];
  o_globalData.vandermondeMatrixInverse = &plasticityGlobalMatrixMem[ seissol::model::plasticity_globalMatrixOffsets[1] ];
}

void seissol::initializers::MemoryManager::allocateIntegrationBufferLTS() {
  /*
   *  (thread-local) LTS integration buffers, allocate
   */
  int l_numberOfThreads = 1;
#ifdef _OPENMP
  l_numberOfThreads = omp_get_max_threads();
#endif
  m_integrationBufferLTS = (real*) m_memoryAllocator.allocateMemory( l_numberOfThreads*(4*NUMBER_OF_ALIGNED_DOFS)*sizeof(real), PAGESIZE_STACK, MEMKIND_GLOBAL ) ;

  /*
   *  (thread-local) LTS integration buffers, initialize
   */
#ifdef _OPENMP
  #pragma omp parallel
  {
    size_t l_threadOffset = omp_get_thread_num()*(4*NUMBER_OF_ALIGNED_DOFS);
#else
    size_t l_threadOffset = 0;
#endif
    for ( unsigned int l_dof = 0; l_dof < (4*NUMBER_OF_ALIGNED_DOFS); l_dof++ ) {
      m_integrationBufferLTS[l_dof + l_threadOffset] = (real)0.0;
    }
#ifdef _OPENMP
  }
#endif
}

void seissol::initializers::MemoryManager::correctGhostRegionSetups()
{
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    Layer& ghost = m_ltsTree.child(tc).child<Ghost>();
    CellLocalInformation* cellInformation = ghost.var(m_lts.cellInformation);
    
    unsigned int l_offset = 0;
    for( unsigned int l_region = 0; l_region < m_meshStructure[tc].numberOfRegions; l_region++ ) {
      // iterate over ghost cells
      for( unsigned int l_cell = 0; l_cell < m_meshStructure[tc].numberOfGhostRegionCells[l_region]; l_cell++ ) {
        if( l_cell < m_meshStructure[tc].numberOfGhostRegionDerivatives[l_region] ) {
          // assert the cell provides derivatives
          assert( (cellInformation[l_offset+l_cell].ltsSetup >> 9)%2 );

          // reset possible buffers
          cellInformation[l_offset+l_cell].ltsSetup &= ( ~(1 << 8 ) );
          cellInformation[l_offset+l_cell].ltsSetup &= ( ~(1 << 10) );
        } else {
          // assert the cell provides buffers
          assert( (cellInformation[l_offset+l_cell].ltsSetup >> 8)%2 );

          // reset possible derivatives
          cellInformation[l_offset+l_cell].ltsSetup &= ( ~(1 << 9 ) );
        }
      }
      // update offset with ghost region size
      l_offset +=  m_meshStructure[tc].numberOfGhostRegionCells[l_region];
    }
  }
}

void seissol::initializers::MemoryManager::deriveLayerLayouts() {
  // initialize memory
#ifdef USE_MPI
  m_numberOfGhostBuffers           = (unsigned int*)  m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int  ), 1 );
  m_numberOfGhostRegionBuffers     = (unsigned int**) m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int* ), 1 );
  m_numberOfGhostDerivatives       = (unsigned int*)  m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int  ), 1 );
  m_numberOfGhostRegionDerivatives = (unsigned int**) m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int* ), 1 );

  m_numberOfCopyBuffers            = (unsigned int*)  m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int  ), 1 );
  m_numberOfCopyRegionBuffers      = (unsigned int**) m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int* ), 1 );
  m_numberOfCopyDerivatives        = (unsigned int*)  m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int  ), 1 );
  m_numberOfCopyRegionDerivatives  = (unsigned int**) m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int* ), 1 );
#endif // USE_MPI

  m_numberOfInteriorBuffers        = (unsigned int*)  m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int  ), 1 );
  m_numberOfInteriorDerivatives    = (unsigned int*)  m_memoryAllocator.allocateMemory( m_ltsTree.numChildren() * sizeof( unsigned int  ), 1 );

  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
#ifdef USE_MPI
    CellLocalInformation* ghostCellInformation    = cluster.child<Ghost>().var(m_lts.cellInformation);
    CellLocalInformation* copyCellInformation     = cluster.child<Copy>().var(m_lts.cellInformation);
#endif
    CellLocalInformation* interiorCellInformation = cluster.child<Interior>().var(m_lts.cellInformation);
#ifdef USE_MPI
    m_numberOfGhostBuffers[             tc] = 0;
    m_numberOfGhostRegionBuffers[       tc] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[tc].numberOfRegions * sizeof( unsigned int ), 1 );
    m_numberOfGhostDerivatives[         tc] = 0;
    m_numberOfGhostRegionDerivatives[   tc] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[tc].numberOfRegions * sizeof( unsigned int ), 1 );

    m_numberOfCopyBuffers[              tc] = 0;
    m_numberOfCopyRegionBuffers[        tc] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[tc].numberOfRegions * sizeof( unsigned int ), 1 );
    m_numberOfCopyDerivatives[          tc] = 0;
    m_numberOfCopyRegionDerivatives[    tc] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[tc].numberOfRegions * sizeof( unsigned int ), 1 );
#endif // USE_MPI

    m_numberOfInteriorBuffers[          tc]       = 0;
    m_numberOfInteriorDerivatives[      tc]       = 0;

#ifdef USE_MPI
    unsigned int l_ghostOffset = 0;
    unsigned int l_copyOffset  = 0;
    for( unsigned int l_region = 0; l_region < m_meshStructure[tc].numberOfRegions; l_region++ ) {
      m_numberOfGhostRegionBuffers[     tc][l_region] = 0;
      m_numberOfGhostRegionDerivatives[ tc][l_region] = 0;

      m_numberOfCopyRegionBuffers[      tc][l_region] = 0;
      m_numberOfCopyRegionDerivatives[  tc][l_region] = 0;

      // iterate over all cells of this clusters ghost layer
      for( unsigned int l_cell = 0; l_cell < m_meshStructure[tc].numberOfGhostRegionCells[l_region]; l_cell++ ) {
        // ensure that either buffers or derivatives are used; not both!
        bool l_buffer      = ( ghostCellInformation[l_cell+l_ghostOffset].ltsSetup >> 8 ) % 2;
        bool l_derivatives = ( ghostCellInformation[l_cell+l_ghostOffset].ltsSetup >> 9 ) % 2;

        if( (l_buffer && l_derivatives) || ( l_buffer || l_derivatives ) == false ) logError() << "invalid ghost lts setup" << l_buffer << l_derivatives;

        // check if this cell requires a buffer and/or derivatives
        if( ( ghostCellInformation[l_cell+l_ghostOffset].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfGhostRegionBuffers[    tc][l_region]++;
        if( ( ghostCellInformation[l_cell+l_ghostOffset].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfGhostRegionDerivatives[tc][l_region]++;
      }
      l_ghostOffset += m_meshStructure[tc].numberOfGhostRegionCells[l_region];

      // iterate over all cells of this clusters copy layer
      for( unsigned int l_cell = 0; l_cell < m_meshStructure[tc].numberOfCopyRegionCells[l_region]; l_cell++ ) {
        // assert that buffers or derivatives are requested
        assert( ( ( copyCellInformation[l_cell+l_copyOffset].ltsSetup >> 8 ) % 2 ||
                  ( copyCellInformation[l_cell+l_copyOffset].ltsSetup >> 9 ) % 2 )
                == true );

        // check if this cell requires a buffer and/or derivatives
        if( ( copyCellInformation[l_cell+l_copyOffset].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfCopyRegionBuffers[    tc][l_region]++;
        if( ( copyCellInformation[l_cell+l_copyOffset].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfCopyRegionDerivatives[tc][l_region]++;
      }
      l_copyOffset += m_meshStructure[tc].numberOfCopyRegionCells[l_region];

      // update number of buffers and derivatives in the ghost and copy layer
      m_numberOfGhostBuffers[     tc ]  += m_numberOfGhostRegionBuffers[     tc][l_region];
      m_numberOfGhostDerivatives[ tc ]  += m_numberOfGhostRegionDerivatives[ tc][l_region];
      m_numberOfCopyBuffers[      tc ]  += m_numberOfCopyRegionBuffers[      tc][l_region];
      m_numberOfCopyDerivatives[  tc ]  += m_numberOfCopyRegionDerivatives[  tc][l_region];
    }
#endif // USE_MPI

    // iterate over all cells of this clusters interior
    for( unsigned int l_cell = 0; l_cell < m_meshStructure[tc].numberOfInteriorCells; l_cell++ ) {
      // check if this cell requires a buffer and/or derivatives
      if( ( interiorCellInformation[l_cell].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfInteriorBuffers[    tc]++;
      if( ( interiorCellInformation[l_cell].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfInteriorDerivatives[tc]++;
    }
  }
}

#ifdef USE_MPI
void seissol::initializers::MemoryManager::initializeCommunicationStructure() {
  // reset mpi requests
  for( unsigned int l_cluster = 0; l_cluster < m_ltsTree.numChildren(); l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      m_meshStructure[l_cluster].sendRequests[l_region] = MPI_REQUEST_NULL;
      m_meshStructure[l_cluster].receiveRequests[l_region] = MPI_REQUEST_NULL;
    }
  }

  /*
   * ghost layer
   */
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
    real* ghostStart = static_cast<real*>(cluster.child<Ghost>().bucket(m_lts.buffersDerivatives));
    for( unsigned int l_region = 0; l_region < m_meshStructure[tc].numberOfRegions; l_region++ ) {
      // set pointer to ghost region
      m_meshStructure[tc].ghostRegions[l_region] = ghostStart;

      // derive the ghost region size
      unsigned int l_numberOfDerivatives = m_meshStructure[tc].numberOfGhostRegionDerivatives[l_region];
      unsigned int l_numberOfBuffers     = m_meshStructure[tc].numberOfGhostRegionCells[l_region] - l_numberOfDerivatives;

      // set size
      m_meshStructure[tc].ghostRegionSizes[l_region] = NUMBER_OF_ALIGNED_DOFS * l_numberOfBuffers +
                                                       NUMBER_OF_ALIGNED_DERS * l_numberOfDerivatives;

      // update the pointer
      ghostStart += m_meshStructure[tc].ghostRegionSizes[l_region];
    }
  }

  /*
   * copy layer
   */   
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    Layer& copy = m_ltsTree.child(tc).child<Copy>();
    real** buffers = copy.var(m_lts.buffers);
    real** derivatives = copy.var(m_lts.derivatives);
    // copy region offset
    unsigned int l_offset = 0;

    for( unsigned int l_region = 0; l_region < m_meshStructure[tc].numberOfRegions; l_region++ ) {
      // derive the communication size
      unsigned int l_numberOfDerivatives = m_meshStructure[tc].numberOfCommunicatedCopyRegionDerivatives[l_region];
      unsigned int l_numberOfBuffers     = m_meshStructure[tc].numberOfCopyRegionCells[l_region] - l_numberOfDerivatives;

      // assert number of communicated buffers fits into total number of buffers
      assert( m_numberOfCopyRegionBuffers[tc][l_region] >= l_numberOfBuffers );

      // set pointer to copy region start
      if( l_numberOfBuffers > 0 ) {
        m_meshStructure[tc].copyRegions[l_region] = buffers[l_numberOfDerivatives + l_offset];
      }
      else {
        m_meshStructure[tc].copyRegions[l_region] = derivatives[l_offset];
      }

      // assert the pointer is set
      assert( m_meshStructure[tc].copyRegions[l_region] != NULL );

      // set size
      m_meshStructure[tc].copyRegionSizes[l_region] = NUMBER_OF_ALIGNED_DOFS * l_numberOfBuffers +
                                                      NUMBER_OF_ALIGNED_DERS * l_numberOfDerivatives;

      // jump over region
      l_offset += m_meshStructure[tc].numberOfCopyRegionCells[l_region];
    }
  }
}
#endif

void seissol::initializers::MemoryManager::initializeFaceNeighbors( unsigned    cluster,
                                                                    Layer&      layer )
{
#ifdef USE_MPI
  assert(layer.getLayerType() == Copy || layer.getLayerType() == Interior);
#else
  assert(layer.getLayerType() == Interior);
#endif
  
  // iterate over clusters
  
  real** buffers = m_ltsTree.var(m_lts.buffers);          // faceNeighborIds are ltsIds and not layer-local
  real** derivatives = m_ltsTree.var(m_lts.derivatives);  // faceNeighborIds are ltsIds and not layer-local
  real *(*faceNeighbors)[4] = layer.var(m_lts.faceNeighbors);
  CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);
    
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    for (unsigned face = 0; face < 4; ++face) {
      if (  cellInformation[cell].faceTypes[face] == regular
         || cellInformation[cell].faceTypes[face] == periodic
         || cellInformation[cell].faceTypes[face] == dynamicRupture ) {
        // neighboring cell provides derivatives
        if( (cellInformation[cell].ltsSetup >> face) % 2 ) {
          faceNeighbors[cell][face] = derivatives[ cellInformation[cell].faceNeighborIds[face] ];
        }
        // neighboring cell provides a time buffer
        else {
          faceNeighbors[cell][face] = buffers[ cellInformation[cell].faceNeighborIds[face] ];
        }
        assert( faceNeighbors[cell][face] != NULL );
      }
      // free surface boundary
      else if( cellInformation[cell].faceTypes[face] == freeSurface ) {
        if( (cellInformation[cell].ltsSetup >> face) % 2 == 0 ) { // free surface on buffers
          faceNeighbors[cell][face] = layer.var(m_lts.buffers)[cell];
        }
        else { // free surface on derivatives
          faceNeighbors[cell][face] = layer.var(m_lts.derivatives)[cell];
        }
        assert( faceNeighbors[cell][face] != NULL );
      }
      // absorbing
      else if( cellInformation[cell].faceTypes[face] == outflow ) {
        // NULL pointer; absorbing: data is not used
        faceNeighbors[cell][face] = NULL;
      }
      else {
        // assert all cases are covered
        assert( false );
      }
    }
  }
}

void seissol::initializers::MemoryManager::initializeBuffersDerivatives() {
  // initialize the pointers of the internal state
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
#ifdef USE_MPI
    /*
     * ghost layer
     */
    InternalState::setUpLayerPointers( m_meshStructure[tc].numberOfRegions,
                                       m_meshStructure[tc].numberOfGhostRegionCells,
                                       cluster.child<Ghost>().var(m_lts.cellInformation),
                                       m_numberOfGhostRegionBuffers[tc],
                                       m_numberOfGhostRegionDerivatives[tc],
                                       static_cast<real*>(cluster.child<Ghost>().bucket(m_lts.buffersDerivatives)),
                                       cluster.child<Ghost>().var(m_lts.buffers),
                                       cluster.child<Ghost>().var(m_lts.derivatives) );

    /*
     * Copy layer
     */
    InternalState::setUpLayerPointers( m_meshStructure[tc].numberOfRegions,
                                       m_meshStructure[tc].numberOfCopyRegionCells,
                                       cluster.child<Copy>().var(m_lts.cellInformation),
                                       m_numberOfCopyRegionBuffers[tc],
                                       m_numberOfCopyRegionDerivatives[tc],
                                       static_cast<real*>(cluster.child<Copy>().bucket(m_lts.buffersDerivatives)),
                                       cluster.child<Copy>().var(m_lts.buffers),
                                       cluster.child<Copy>().var(m_lts.derivatives) );
#endif

    /*
     * Interior
     */
    InternalState::setUpInteriorPointers( m_meshStructure[tc].numberOfInteriorCells,
                                          cluster.child<Interior>().var(m_lts.cellInformation),
                                          m_numberOfInteriorBuffers[tc],
                                          m_numberOfInteriorDerivatives[tc],
                                          static_cast<real*>(cluster.child<Interior>().bucket(m_lts.buffersDerivatives)),
                                          cluster.child<Interior>().var(m_lts.buffers),
                                          cluster.child<Interior>().var(m_lts.derivatives)  );
  }
}

void seissol::initializers::MemoryManager::touchBuffersDerivatives( Layer& layer ) {
  real** buffers = layer.var(m_lts.buffers);
  real** derivatives = layer.var(m_lts.derivatives);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    // touch buffers
    real* buffer = buffers[cell];
    if (buffer != NULL) {
      for (unsigned dof = 0; dof < NUMBER_OF_ALIGNED_DOFS; ++dof) {
          // zero time integration buffers
          buffer[dof] = (real) 0;
      }
    }

    // touch derivatives
    real* derivative = derivatives[cell];
    if (derivative != NULL) {
      for (unsigned dof = 0; dof < NUMBER_OF_ALIGNED_DERS; ++dof ) {
        derivative[dof] = (real) 0;
      }
    }
  }
}

void seissol::initializers::MemoryManager::fixateLtsTree( struct TimeStepping&        i_timeStepping,
                                                          struct MeshStructure*       i_meshStructure,
                                                          unsigned                    globalDynamicRuptureTimeCluster,
                                                          unsigned                    numberOfDRCopyFaces,
                                                          unsigned                    numberOfDRInteriorFaces )
{
  // store mesh structure and the number of time clusters
  m_meshStructure = i_meshStructure;
  
  // Setup tree variables
  m_lts.addTo(m_ltsTree);
  seissol::SeisSol::main.postProcessor().allocateMemory(&m_ltsTree);
  m_ltsTree.setNumberOfTimeClusters(i_timeStepping.numberOfLocalClusters);

  /// From this point, the tree layout, variables, and buckets cannot be changed anymore
  m_ltsTree.fixate();
  
  // Set number of cells and bucket sizes in ltstree
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
    cluster.child<Ghost>().setNumberOfCells(i_meshStructure[tc].numberOfGhostCells);
    cluster.child<Copy>().setNumberOfCells(i_meshStructure[tc].numberOfCopyCells);
    cluster.child<Interior>().setNumberOfCells(i_meshStructure[tc].numberOfInteriorCells);
  }

  m_ltsTree.allocateVariables();
  m_ltsTree.touchVariables();
  
  /// Dynamic rupture tree  
  m_dynRup.addTo(m_dynRupTree);
  m_dynRupTree.setNumberOfTimeClusters(i_timeStepping.numberOfLocalClusters);
  m_dynRupTree.fixate();

  unsigned localDynamicRuptureTimeCluster = std::numeric_limits<unsigned>::max();
  assert( i_timeStepping.numberOfLocalClusters < localDynamicRuptureTimeCluster);
  for (unsigned cluster = 0; cluster < i_timeStepping.numberOfLocalClusters; ++cluster) {
    if (i_timeStepping.clusterIds[cluster] == globalDynamicRuptureTimeCluster) {
      localDynamicRuptureTimeCluster = cluster;
    }
  }

  for (unsigned tc = 0; tc < m_dynRupTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_dynRupTree.child(tc);
    cluster.child<Ghost>().setNumberOfCells(0);
    if (tc == localDynamicRuptureTimeCluster) {
      cluster.child<Copy>().setNumberOfCells(numberOfDRCopyFaces);
      cluster.child<Interior>().setNumberOfCells(numberOfDRInteriorFaces);
    } else {
      cluster.child<Copy>().setNumberOfCells(0);
      cluster.child<Interior>().setNumberOfCells(0);
    }
  }
  
  m_dynRupTree.allocateVariables();
  m_dynRupTree.touchVariables();
}

void seissol::initializers::MemoryManager::deriveDisplacementsBucket()
{
  for ( seissol::initializers::LTSTree::leaf_iterator layer = m_ltsTree.beginLeaf(m_lts.displacements.mask); layer != m_ltsTree.endLeaf(); ++layer) {
    CellLocalInformation* cellInformation = layer->var(m_lts.cellInformation);
    real** displacements = layer->var(m_lts.displacements);
    
    unsigned numberOfCells = 0;
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      bool hasFreeSurface = false;
      for (unsigned face = 0; face < 4; ++face) {
        hasFreeSurface = hasFreeSurface || (cellInformation[cell].faceTypes[face] == freeSurface);
      }
      if (hasFreeSurface) {
        // We add the base address later when the bucket is allocated
        // +1 is necessary as we want to reserve the NULL pointer for cell without displacement.
        displacements[cell] = static_cast<real*>(NULL) + 1 + numberOfCells * NUMBER_OF_ALIGNED_VELOCITY_DOFS;
        ++numberOfCells;
      } else {
        displacements[cell] = NULL;
      }
    }
    layer->setBucketSize(m_lts.displacementsBuffer, numberOfCells * NUMBER_OF_ALIGNED_VELOCITY_DOFS * sizeof(real));
  }
}

void seissol::initializers::MemoryManager::initializeDisplacements()
{
  for ( seissol::initializers::LTSTree::leaf_iterator layer = m_ltsTree.beginLeaf(m_lts.displacements.mask); layer != m_ltsTree.endLeaf(); ++layer) {
    real** displacements = layer->var(m_lts.displacements);
    real* bucket = static_cast<real*>(layer->bucket(m_lts.displacementsBuffer));

    #pragma omp parallel for schedule(static)
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      if (displacements[cell] != NULL) {
        displacements[cell] = bucket + ((displacements[cell] - static_cast<real*>(NULL)) - 1);
        for (unsigned dof = 0; dof < NUMBER_OF_ALIGNED_VELOCITY_DOFS; ++dof) {
          // zero displacements
          displacements[cell][dof] = static_cast<real>(0.0);
        }
      }
    }
  }
}

void seissol::initializers::MemoryManager::initializeMemoryLayout(bool enableFreeSurfaceIntegration)
{
  // correct LTS-information in the ghost layer
  correctGhostRegionSetups();

  // derive the layouts of the layers
  deriveLayerLayouts();

  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);

    size_t l_ghostSize = 0;
    size_t l_copySize = 0;
    size_t l_interiorSize = 0;
#ifdef USE_MPI
    for( unsigned int l_region = 0; l_region < m_meshStructure[tc].numberOfRegions; l_region++ ) {
      l_ghostSize    += sizeof(real) * NUMBER_OF_ALIGNED_DOFS * m_numberOfGhostRegionBuffers[tc][l_region];
      l_ghostSize    += sizeof(real) * NUMBER_OF_ALIGNED_DERS * m_numberOfGhostRegionDerivatives[tc][l_region];

      l_copySize     += sizeof(real) * NUMBER_OF_ALIGNED_DOFS * m_numberOfCopyRegionBuffers[tc][l_region];
      l_copySize     += sizeof(real) * NUMBER_OF_ALIGNED_DERS * m_numberOfCopyRegionDerivatives[tc][l_region];
    }
#endif // USE_MPI
    l_interiorSize += sizeof(real) * NUMBER_OF_ALIGNED_DOFS * m_numberOfInteriorBuffers[tc];
    l_interiorSize += sizeof(real) * NUMBER_OF_ALIGNED_DERS * m_numberOfInteriorDerivatives[tc];
    
    cluster.child<Ghost>().setBucketSize(m_lts.buffersDerivatives, l_ghostSize);
    cluster.child<Copy>().setBucketSize(m_lts.buffersDerivatives, l_copySize);
    cluster.child<Interior>().setBucketSize(m_lts.buffersDerivatives, l_interiorSize);
  }
  
  if (enableFreeSurfaceIntegration) {
    deriveDisplacementsBucket();
  }
  
  
  m_ltsTree.allocateBuckets();

  // initialize the internal state
  initializeBuffersDerivatives();

  // initialize face neighbors
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
#ifdef USE_MPI
    initializeFaceNeighbors(tc, cluster.child<Copy>());
#endif
    initializeFaceNeighbors(tc, cluster.child<Interior>());
  }

  for (LTSTree::leaf_iterator it = m_ltsTree.beginLeaf(); it != m_ltsTree.endLeaf(); ++it) {
    touchBuffersDerivatives(*it);
  }

#ifdef USE_MPI
  // initialize the communication structure
  initializeCommunicationStructure();
#endif
  
  if (enableFreeSurfaceIntegration) {
    initializeDisplacements();
  }
}

void seissol::initializers::MemoryManager::getMemoryLayout( unsigned int                    i_cluster,
                                                            struct MeshStructure          *&o_meshStructure,
                                                            struct GlobalData             *&o_globalData
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
                                                            ,struct GlobalData             *&o_globalDataCopies
#endif
                                                          ) {
  o_meshStructure           =  m_meshStructure + i_cluster;
  o_globalData              = &m_globalData;
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
  o_globalDataCopies        =  m_globalDataCopies;
#endif
}
