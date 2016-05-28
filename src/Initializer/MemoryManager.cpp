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
#include "MemoryManager.h"
#include "InternalState.h"

#include <Kernels/common.hpp>

// if equations == viscoelastic
#ifdef REQUIRE_SOURCE_MATRIX
#include <generated_code/init.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if CONVERGENCE_ORDER == 2
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#define MEMKIND_CONSTANT seissol::memory::HighBandwidth
#define MEMKIND_DOFS     seissol::memory::HighBandwidth
#elif CONVERGENCE_ORDER == 3
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#define MEMKIND_CONSTANT seissol::memory::HighBandwidth
#define MEMKIND_DOFS     seissol::memory::HighBandwidth
#elif CONVERGENCE_ORDER == 4
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#define MEMKIND_CONSTANT seissol::memory::HighBandwidth
#define MEMKIND_DOFS     seissol::memory::Standard
#elif CONVERGENCE_ORDER == 5
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#define MEMKIND_CONSTANT seissol::memory::Standard
#define MEMKIND_DOFS     seissol::memory::Standard
#elif CONVERGENCE_ORDER == 6
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#define MEMKIND_CONSTANT seissol::memory::Standard
#define MEMKIND_DOFS     seissol::memory::Standard
#elif CONVERGENCE_ORDER == 7
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#define MEMKIND_CONSTANT seissol::memory::Standard
#define MEMKIND_DOFS     seissol::memory::Standard
#elif CONVERGENCE_ORDER == 8
#define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#define MEMKIND_TIMEDOFS seissol::memory::Standard
#define MEMKIND_CONSTANT seissol::memory::Standard
#define MEMKIND_DOFS     seissol::memory::Standard
#else
#error Preprocessor flag CONVERGENCE_ORDER is not in {2, 3, 4, 5, 6, 7, 8}.
#endif

seissol::initializers::MemoryManager::MemoryManager()
  : m_integrationBufferLTS(NULL)
{
}

void seissol::initializers::MemoryManager::initialize( const seissol::XmlParser &i_matrixReader )
{
#ifndef REQUIRE_SOURCE_MATRIX
  // init the sparse switch
#   define SPARSE_SWITCH
#   include <initialization/bind.h>
#   undef SPARSE_SWITCH
#endif

  // allocate thread-local LTS integration buffers
  allocateIntegrationBufferLTS();
  
// if equations == viscoelastic
// @TODO Remove ifdef and generalize initialization
// @TODO This implementation doesn't backport the support for multiple copies of global data
// @TODO make sure that multiple copies of global data are still supported by this unified implementation
#ifdef REQUIRE_SOURCE_MATRIX
  real* globalMatrixMem = static_cast<real*>(m_memoryAllocator.allocateMemory( seissol::model::globalMatrixOffsets[seissol::model::numGlobalMatrices] * sizeof(real), PAGESIZE_HEAP, MEMKIND_GLOBAL ));
  for (unsigned matrix = 0; matrix < seissol::model::numGlobalMatrices; ++matrix) {
    memcpy(
      &globalMatrixMem[ seissol::model::globalMatrixOffsets[matrix] ],
      seissol::model::globalMatrixValues[matrix],
      (seissol::model::globalMatrixOffsets[matrix+1] - seissol::model::globalMatrixOffsets[matrix]) * sizeof(real)
    );
  }
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    m_globalData.stiffnessMatricesTransposed[transposedStiffness] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[transposedStiffness] ];
  }
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    m_globalData.stiffnessMatrices[stiffness] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[3 + stiffness] ];
  }
  for (unsigned flux = 0; flux < 52; ++flux) {
    m_globalData.fluxMatrices[flux] = &globalMatrixMem[ seissol::model::globalMatrixOffsets[6 + flux] ];
  }

  // @TODO Integrate this step into the code generator
  for (unsigned transposedStiffness = 0; transposedStiffness < 3; ++transposedStiffness) {
    real* matrix = &globalMatrixMem[ seissol::model::globalMatrixOffsets[transposedStiffness] ];
    for (unsigned i = 0; i < seissol::model::globalMatrixOffsets[transposedStiffness+1]-seissol::model::globalMatrixOffsets[transposedStiffness]; ++i) {
      matrix[i] *= -1.0;
    }
  }

  // set LTS integration buffer
  m_globalData.integrationBufferLTS = m_integrationBufferLTS;
#else
  // initialize global matrices
#ifndef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
  initializeGlobalMatrices( i_matrixReader, m_globalData );
#else
  // determine the number of threads
  unsigned int l_numberOfThreads = omp_get_max_threads();
  unsigned int l_numberOfCopiesCeil = (l_numberOfThreads%NUMBER_OF_THREADS_PER_GLOBALDATA_COPY == 0) ? 0 : 1;
  unsigned int l_numberOfCopies = (l_numberOfThreads/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY) + l_numberOfCopiesCeil;
  logInfo(0) << "Number of GlobalData copies: " << l_numberOfCopies;

  m_globalDataCopies = new GlobalData[l_numberOfCopies]; 

  // initialize all copies
#if 0 
  // use serial initialization -> first touch places everything in the corresponding NUMA nodes
  for ( unsigned int l_globalDataCopy = 0; l_globalDataCopy < l_numberOfCopies; l_globalDataCopy++ ) {
    initializeGlobalMatrices( i_matrixReader, m_globalDataCopies[l_globalDataCopy] );
  }
#else
  // initialize in parallel to obtain best possible NUMA placement
  #pragma omp parallel
  {
    if (omp_get_thread_num()%NUMBER_OF_THREADS_PER_GLOBALDATA_COPY == 0) {
      // @TODO check why initializeGlobalMatrices is not thread-safe
      #pragma omp critical
      {
        initializeGlobalMatrices( i_matrixReader, m_globalDataCopies[omp_get_thread_num()/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY] );
      }
    }
  }
#endif
  // set master structure
  m_globalData = m_globalDataCopies[0];
#endif
#endif
}

seissol::initializers::MemoryManager::~MemoryManager() {
  // free members
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
  delete[] m_globalDataCopies;
#endif
}

#ifndef REQUIRE_SOURCE_MATRIX
void seissol::initializers::MemoryManager::initializeGlobalMatrix(          int                        i_sparse,
                                                                   unsigned int                        i_leadingDimension,
                                                                   unsigned int                        i_numberOfColumns,
                                                                      const std::vector<unsigned int> &i_rows,
                                                                      const std::vector<unsigned int> &i_columns,
                                                                      const std::vector<double>       &i_values,
                                                                            real*                      o_matrix ) {
  // assert matching dimensions
  assert( i_rows.size()    == i_columns.size() );
  assert( i_columns.size() == i_values.size()  );

  // easy case: just write the values one after another
  if( i_sparse != -1 ) {
    // assert we have all nonzeros
    assert( static_cast<int>(i_values.size()) == i_sparse );

    for( unsigned int l_entry = 0; l_entry < i_values.size(); l_entry++) {
      o_matrix[l_entry] = i_values[l_entry];
    }
  }
  // dense matrix: set everything to zero and set only nonzeros
  else {
    // set everything to zero
    std::fill( o_matrix, o_matrix+i_leadingDimension*i_numberOfColumns, 0 );

    // iterate over nonzeros
    for( unsigned int l_entry = 0; l_entry < i_values.size(); l_entry++) {
      // index calculation (counting in XML starts at 1)
      unsigned int l_row    = i_rows[l_entry]    - 1;
      unsigned int l_column = i_columns[l_entry] - 1;

      // assert a valid position in the (possibly reduced) size of the matrix
      assert( l_row < i_leadingDimension );
      assert( l_column < i_numberOfColumns ) ;

      // jump over columns
      unsigned int l_position = l_column * i_leadingDimension;
      // jump over rows
      l_position += l_row; 

      // set the nonzero
      o_matrix[l_position] = i_values[l_entry];
    } 
  }
}

void seissol::initializers::MemoryManager::initializeGlobalMatrices( const seissol::XmlParser &i_matrixReader,
                                                                     struct GlobalData        &o_globalData ) {
  /*
   * Test whether LTS integation buffer was allocated
   **/
  if ( m_integrationBufferLTS == NULL ) {
    logError() << "MemoryManager: allocateIntegrationBufferLTS need called before initializeGlobalMatrices!"; 
  }

  /*
   * read the global matrices
   */
  //! vectors, which hold information about our matrices
  std::vector< unsigned int > l_matrixIds;
  std::vector< std::string  > l_matrixNames;
  std::vector< unsigned int > l_matrixNumberOfRows;
  std::vector< unsigned int > l_matrixNumberOfColumns;
  std::vector< bool         > l_matrixSparsities;

  // element information in coordinate format
  std::vector< std::vector<unsigned int> > l_matrixRows;
  std::vector< std::vector<unsigned int> > l_matrixColumns;
  std::vector< std::vector<double>       > l_matrixValues;

  // read the flux matrices
  i_matrixReader.readGlobalMatrices( "flux",
                                     l_matrixIds,  l_matrixNames,
                                     l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                     l_matrixRows, l_matrixColumns, l_matrixValues );
  // assert we have all flux matrices
  assert( l_matrixIds.size() == 52 );

  // read the stiffness matrices
  i_matrixReader.readGlobalMatrices( "stiffness",
                                     l_matrixIds,  l_matrixNames,
                                     l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                     l_matrixRows, l_matrixColumns, l_matrixValues );

  // assert we have all stiffness matrices
  assert( l_matrixIds.size() == 58 );

  // set negative sign of stiffness matrices in time derivative computation
  for( unsigned int l_transposedStiff = 55; l_transposedStiff < 58; l_transposedStiff++ ) {
    for( unsigned int l_entry = 0; l_entry < l_matrixValues[l_transposedStiff].size(); l_entry++ ) {
      l_matrixValues[l_transposedStiff][l_entry] = -l_matrixValues[l_transposedStiff][l_entry];
    }
  }
  
  // read the inverse mass matrix
  i_matrixReader.readGlobalMatrices( "inverseMass",
                                     l_matrixIds, l_matrixNames,
                                     l_matrixNumberOfRows, l_matrixNumberOfColumns, l_matrixSparsities,
                                     l_matrixRows, l_matrixColumns, l_matrixValues );

  // assert we have the mass matrix
  assert ( l_matrixIds.size() == 59 );

  /*
   * Allocate memory.
   */
  // offset to reach each of the matrices; ordering transposed stiffness, stiffness, flux, inverse mass, dummy for total size
  unsigned int l_offset[60]; l_offset[0] = 0;

  // number of aligned reals of a matrix
  unsigned int l_alignedReals = 0;

  // transposed stiffness matrices
  for( unsigned int l_matrix = 0; l_matrix < 3; l_matrix++ ) {
    // dense
    if( m_sparseSwitch[l_matrix+56] == -1 ) {
      l_alignedReals = seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-1 ) * NUMBER_OF_BASIS_FUNCTIONS;
    }
    // sparse
    else {
      l_alignedReals = seissol::kernels::getNumberOfAlignedReals( m_sparseSwitch[l_matrix+56] );
    }

    l_offset[l_matrix+1] = l_offset[l_matrix+0] + l_alignedReals;
  }

  // stiffness matrices
  for( unsigned int l_matrix = 0; l_matrix < 3; l_matrix++ ) {
    // dense
    if( m_sparseSwitch[l_matrix+53] == -1 ) {
      l_alignedReals = NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * seissol::kernels::getNumberOfBasisFunctions( CONVERGENCE_ORDER-1 );
    }
    // sparse
    else {
      l_alignedReals = seissol::kernels::getNumberOfAlignedReals( m_sparseSwitch[l_matrix+53] );
    }

    l_offset[l_matrix+4] = l_offset[l_matrix+3] + l_alignedReals;
  }

  // flux matrices
  for( unsigned int l_matrix = 0; l_matrix < 52; l_matrix++ ) {
    // dense
    if( m_sparseSwitch[l_matrix] == -1 ) {
      l_alignedReals = NUMBER_OF_ALIGNED_BASIS_FUNCTIONS * NUMBER_OF_BASIS_FUNCTIONS;
    }
    // sparse
    else {
      l_alignedReals = seissol::kernels::getNumberOfAlignedReals( m_sparseSwitch[l_matrix] );
    }

    l_offset[l_matrix+7] = l_offset[l_matrix+6] + l_alignedReals;
  }
  
  // inverse mass matrix
  l_alignedReals = seissol::kernels::getNumberOfAlignedReals( NUMBER_OF_BASIS_FUNCTIONS );
  l_offset[59] = l_offset[58] + l_alignedReals;

  real* l_pointer = (real*) m_memoryAllocator.allocateMemory( l_offset[59] * sizeof(real), PAGESIZE_HEAP, MEMKIND_GLOBAL );

  // init data fast to get contiguous physical memory, or at least increase the chances
  for( unsigned int l_init = 0; l_init < l_offset[59]; l_init++) {
    l_pointer[l_init] = static_cast<real>(0.0);
  }

  /*
   * Set up pointers.
   */
  for( unsigned int l_transposedStiffnessMatrix = 0; l_transposedStiffnessMatrix < 3; l_transposedStiffnessMatrix++ ) {
    o_globalData.stiffnessMatricesTransposed[l_transposedStiffnessMatrix] = l_pointer + l_offset[l_transposedStiffnessMatrix];
  }

  for( unsigned int l_stiffnessMatrix = 0; l_stiffnessMatrix < 3; l_stiffnessMatrix++ ) {
    o_globalData.stiffnessMatrices[l_stiffnessMatrix] = l_pointer + l_offset[l_stiffnessMatrix+3];
  }

  for( unsigned int l_fluxMatrix = 0; l_fluxMatrix < 52; l_fluxMatrix++ ) {
    o_globalData.fluxMatrices[l_fluxMatrix] = l_pointer + l_offset[l_fluxMatrix+6];
  }

  o_globalData.inverseMassMatrix = l_pointer + l_offset[58];

  /**
   * Initialize transposed stiffness matrices.
   **/
  for( int l_transposedStiffnessMatrix = 0; l_transposedStiffnessMatrix < 3; l_transposedStiffnessMatrix++ ) {
    // jump over the 52 flux matrices, flu solver and 3 stiffness matrices
    int l_globalMatrix = l_transposedStiffnessMatrix + 56;

    // initialize the stiffness matrix
    initializeGlobalMatrix( m_sparseSwitch[l_globalMatrix],
                            seissol::kernels::getNumberOfAlignedBasisFunctions( CONVERGENCE_ORDER-1 ),
                            NUMBER_OF_BASIS_FUNCTIONS,
                            l_matrixRows[l_globalMatrix-1], // -1: flux solver is not part of the matrices read from XML
                            l_matrixColumns[l_globalMatrix-1],
                            l_matrixValues[l_globalMatrix-1],
                            o_globalData.stiffnessMatricesTransposed[l_transposedStiffnessMatrix] );
  }

  /*
   * Initialize stiffness matrices.
   */
  for( int l_stiffnessMatrix = 0; l_stiffnessMatrix < 3; l_stiffnessMatrix++ ) {
    // jump over the 52 flux matrices and flux solver
    int l_globalMatrix = l_stiffnessMatrix + 53;

    // initialize the stiffness matrix
    initializeGlobalMatrix( m_sparseSwitch[l_globalMatrix],
                            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                            seissol::kernels::getNumberOfBasisFunctions( CONVERGENCE_ORDER-1 ),
                            l_matrixRows[l_globalMatrix-1],
                            l_matrixColumns[l_globalMatrix-1],
                            l_matrixValues[l_globalMatrix-1],
                            o_globalData.stiffnessMatrices[l_stiffnessMatrix] );
  }

  /*
   * Initialize flux matrices.
   */
  for( int l_fluxMatrix = 0; l_fluxMatrix < 52; l_fluxMatrix++) {
    // initialize the stiffness matrix
    initializeGlobalMatrix( m_sparseSwitch[l_fluxMatrix],
                            NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                            NUMBER_OF_BASIS_FUNCTIONS,
                            l_matrixRows[l_fluxMatrix],
                            l_matrixColumns[l_fluxMatrix],
                            l_matrixValues[l_fluxMatrix],
                            o_globalData.fluxMatrices[l_fluxMatrix] );
  }
  
  /*
   * Initialize inverse mass matrix.
   */
  initializeGlobalMatrix( NUMBER_OF_BASIS_FUNCTIONS,
                          NUMBER_OF_ALIGNED_BASIS_FUNCTIONS,
                          NUMBER_OF_BASIS_FUNCTIONS,
                          l_matrixRows[58],
                          l_matrixColumns[58],
                          l_matrixValues[58],
                          o_globalData.inverseMassMatrix );

  /*
   *  (thread-local) LTS integration buffers
   */
  o_globalData.integrationBufferLTS = m_integrationBufferLTS;
}
#endif

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

void seissol::initializers::MemoryManager::setUpLayers( struct CellLocalInformation *i_cellLocalInformation ) {
  // set total numbers
  m_totalNumberOfGhostCells    = 0;
  m_totalNumberOfCopyCells     = 0;
  m_totalNumberOfInteriorCells = 0;

  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    m_totalNumberOfGhostCells    += m_meshStructure[l_cluster].numberOfGhostCells;
    m_totalNumberOfCopyCells     += m_meshStructure[l_cluster].numberOfCopyCells;
    m_totalNumberOfInteriorCells += m_meshStructure[l_cluster].numberOfInteriorCells;
  }

  m_totalNumberOfCells = m_totalNumberOfGhostCells + m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells;

  // set up cell information pointers
#ifdef USE_MPI
  m_ghostCellInformation    = (CellLocalInformation**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( CellLocalInformation*), 1 );
  m_copyCellInformation     = (CellLocalInformation**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( CellLocalInformation*), 1 );
#endif // USE_MPI
  m_interiorCellInformation = (CellLocalInformation**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( CellLocalInformation*), 1 );

  unsigned int l_offset = 0;
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
#ifdef USE_MPI
    // set up ghost
    m_ghostCellInformation[l_cluster] = i_cellLocalInformation + l_offset;
    l_offset += m_meshStructure[l_cluster].numberOfGhostCells;

    // set up copy
    m_copyCellInformation[l_cluster]  = i_cellLocalInformation + l_offset;
    l_offset += m_meshStructure[l_cluster].numberOfCopyCells;
#endif // USE_MPI

    // set up interior
    m_interiorCellInformation[l_cluster] = i_cellLocalInformation + l_offset;
    l_offset += m_meshStructure[l_cluster].numberOfInteriorCells;
  }
}

void seissol::initializers::MemoryManager::correctGhostRegionSetups( struct CellLocalInformation *io_cellLocalInformation ) {
  unsigned int l_offset = 0;

  // iterate over time clusters
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    // iterate over ghost region
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      // iterate over ghost cells
      for( unsigned int l_cell = 0; l_cell < m_meshStructure[l_cluster].numberOfGhostRegionCells[l_region]; l_cell++ ) {
        if( l_cell < m_meshStructure[l_cluster].numberOfGhostRegionDerivatives[l_region] ) {
          // assert the cell provides derivatives
          assert( (io_cellLocalInformation[l_offset+l_cell].ltsSetup >> 9)%2 );

          // reset possible buffers
          io_cellLocalInformation[l_offset+l_cell].ltsSetup &= ( ~(1 << 8 ) );
          io_cellLocalInformation[l_offset+l_cell].ltsSetup &= ( ~(1 << 10) );
        }
        else {
          // assert the cell provides buffers
          assert( (io_cellLocalInformation[l_offset+l_cell].ltsSetup >> 8)%2 );

          // reset possible derivatives
          io_cellLocalInformation[l_offset+l_cell].ltsSetup &= ( ~(1 << 9 ) );
        }
      }
      // update offset with ghost region size
      l_offset +=  m_meshStructure[l_cluster].numberOfGhostRegionCells[l_region];
    }

    // jump over copy and interior cells
    l_offset += m_meshStructure[l_cluster].numberOfCopyCells;
    l_offset += m_meshStructure[l_cluster].numberOfInteriorCells;
  }
}

void seissol::initializers::MemoryManager::deriveLayerLayouts() {
  // initialize memory
#ifdef USE_MPI
  m_numberOfGhostBuffers           = (unsigned int*)  m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int  ), 1 );
  m_numberOfGhostRegionBuffers     = (unsigned int**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int* ), 1 );
  m_numberOfGhostDerivatives       = (unsigned int*)  m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int  ), 1 );
  m_numberOfGhostRegionDerivatives = (unsigned int**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int* ), 1 );

  m_numberOfCopyBuffers            = (unsigned int*)  m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int  ), 1 );
  m_numberOfCopyRegionBuffers      = (unsigned int**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int* ), 1 );
  m_numberOfCopyDerivatives        = (unsigned int*)  m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int  ), 1 );
  m_numberOfCopyRegionDerivatives  = (unsigned int**) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int* ), 1 );
#endif // USE_MPI

  m_numberOfInteriorBuffers        = (unsigned int*)  m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int  ), 1 );
  m_numberOfInteriorDerivatives    = (unsigned int*)  m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( unsigned int  ), 1 );

  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
#ifdef USE_MPI
    m_numberOfGhostBuffers[             l_cluster] = 0;
    m_numberOfGhostRegionBuffers[       l_cluster] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[l_cluster].numberOfRegions * sizeof( unsigned int ), 1 );
    m_numberOfGhostDerivatives[         l_cluster] = 0;
    m_numberOfGhostRegionDerivatives[   l_cluster] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[l_cluster].numberOfRegions * sizeof( unsigned int ), 1 );

    m_numberOfCopyBuffers[              l_cluster] = 0;
    m_numberOfCopyRegionBuffers[        l_cluster] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[l_cluster].numberOfRegions * sizeof( unsigned int ), 1 );
    m_numberOfCopyDerivatives[          l_cluster] = 0;
    m_numberOfCopyRegionDerivatives[    l_cluster] = (unsigned int*)  m_memoryAllocator.allocateMemory( m_meshStructure[l_cluster].numberOfRegions * sizeof( unsigned int ), 1 );
#endif // USE_MPI

    m_numberOfInteriorBuffers[          l_cluster]       = 0;
    m_numberOfInteriorDerivatives[      l_cluster]       = 0;

#ifdef USE_MPI
    unsigned int l_ghostOffset = 0;
    unsigned int l_copyOffset  = 0;
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      m_numberOfGhostRegionBuffers[     l_cluster][l_region] = 0;
      m_numberOfGhostRegionDerivatives[ l_cluster][l_region] = 0;

      m_numberOfCopyRegionBuffers[      l_cluster][l_region] = 0;
      m_numberOfCopyRegionDerivatives[  l_cluster][l_region] = 0;

      // iterate over all cells of this clusters ghost layer
      for( unsigned int l_cell = 0; l_cell < m_meshStructure[l_cluster].numberOfGhostRegionCells[l_region]; l_cell++ ) {
        // ensure that either buffers or derivatives are used; not both!
        bool l_buffer      = ( m_ghostCellInformation[l_cluster][l_cell+l_ghostOffset].ltsSetup >> 8 ) % 2;
        bool l_derivatives = ( m_ghostCellInformation[l_cluster][l_cell+l_ghostOffset].ltsSetup >> 9 ) % 2;

        if( (l_buffer && l_derivatives) || ( l_buffer || l_derivatives ) == false ) logError() << "invalid ghost lts setup" << l_buffer << l_derivatives;

        // check if this cell requires a buffer and/or derivatives
        if( ( m_ghostCellInformation[l_cluster][l_cell+l_ghostOffset].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfGhostRegionBuffers[    l_cluster][l_region]++;
        if( ( m_ghostCellInformation[l_cluster][l_cell+l_ghostOffset].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfGhostRegionDerivatives[l_cluster][l_region]++;
      }
      l_ghostOffset += m_meshStructure[l_cluster].numberOfGhostRegionCells[l_region];

      // iterate over all cells of this clusters copy layer
      for( unsigned int l_cell = 0; l_cell < m_meshStructure[l_cluster].numberOfCopyRegionCells[l_region]; l_cell++ ) {
        // assert that buffers or derivatives are requested
        assert( ( ( m_copyCellInformation[l_cluster][l_cell+l_copyOffset].ltsSetup >> 8 ) % 2 ||
                  ( m_copyCellInformation[l_cluster][l_cell+l_copyOffset].ltsSetup >> 9 ) % 2 )
                == true );

        // check if this cell requires a buffer and/or derivatives
        if( ( m_copyCellInformation[l_cluster][l_cell+l_copyOffset].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfCopyRegionBuffers[    l_cluster][l_region]++;
        if( ( m_copyCellInformation[l_cluster][l_cell+l_copyOffset].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfCopyRegionDerivatives[l_cluster][l_region]++;
      }
      l_copyOffset += m_meshStructure[l_cluster].numberOfCopyRegionCells[l_region];

      // update number of buffers and derivatives in the ghost and copy layer
      m_numberOfGhostBuffers[     l_cluster ]  += m_numberOfGhostRegionBuffers[     l_cluster][l_region];
      m_numberOfGhostDerivatives[ l_cluster ]  += m_numberOfGhostRegionDerivatives[ l_cluster][l_region];
      m_numberOfCopyBuffers[      l_cluster ]  += m_numberOfCopyRegionBuffers[      l_cluster][l_region];
      m_numberOfCopyDerivatives[  l_cluster ]  += m_numberOfCopyRegionDerivatives[  l_cluster][l_region];
    }
#endif // USE_MPI

    // iterate over all cells of this clusters interior
    for( unsigned int l_cell = 0; l_cell < m_meshStructure[l_cluster].numberOfInteriorCells; l_cell++ ) {
      // check if this cell requires a buffer and/or derivatives
      if( ( m_interiorCellInformation[l_cluster][l_cell].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfInteriorBuffers[    l_cluster]++;
      if( ( m_interiorCellInformation[l_cluster][l_cell].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfInteriorDerivatives[l_cluster]++;
    }
  }
}

#ifdef USE_MPI
void seissol::initializers::MemoryManager::initializeCommunicationStructure() {
  // reset mpi requests
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      m_meshStructure[l_cluster].sendRequests[l_region] = MPI_REQUEST_NULL;
      m_meshStructure[l_cluster].receiveRequests[l_region] = MPI_REQUEST_NULL;
    }
  }

  /*
   * ghost layer
   */
  real *l_ghostStart = m_internalState.ghostLayer;
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      // set pointer to ghost region
      m_meshStructure[l_cluster].ghostRegions[l_region] = l_ghostStart;

      // derive the ghost region size
      unsigned int l_numberOfDerivatives = m_meshStructure[l_cluster].numberOfGhostRegionDerivatives[l_region];
      unsigned int l_numberOfBuffers     = m_meshStructure[l_cluster].numberOfGhostRegionCells[l_region] - l_numberOfDerivatives;

      // set size
      m_meshStructure[l_cluster].ghostRegionSizes[l_region] = NUMBER_OF_ALIGNED_DOFS * l_numberOfBuffers +
                                                              NUMBER_OF_ALIGNED_DERS * l_numberOfDerivatives;

      // update the pointer
      l_ghostStart += m_meshStructure[l_cluster].ghostRegionSizes[l_region];
    }
  }

  /*
   * copy layer
   */
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    // copy region offset
    unsigned int l_offset = 0;

    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      // derive the communication size
      unsigned int l_numberOfDerivatives = m_meshStructure[l_cluster].numberOfCommunicatedCopyRegionDerivatives[l_region];
      unsigned int l_numberOfBuffers     = m_meshStructure[l_cluster].numberOfCopyRegionCells[l_region] - l_numberOfDerivatives;

      // assert number of communicated buffers fits into total number of buffers
      assert( m_numberOfCopyRegionBuffers[l_cluster][l_region] >= l_numberOfBuffers );

      // set pointer to copy region start
      if( l_numberOfBuffers > 0 ) {
        m_meshStructure[l_cluster].copyRegions[l_region] = m_cells[l_cluster].copyBuffers[l_numberOfDerivatives+l_offset];
      }
      else {
        m_meshStructure[l_cluster].copyRegions[l_region] = m_cells[l_cluster].copyDerivatives[l_offset];
      }

      // assert the pointer is set
      assert( m_meshStructure[l_cluster].copyRegions[l_region] != NULL );

      // set size
      m_meshStructure[l_cluster].copyRegionSizes[l_region] = NUMBER_OF_ALIGNED_DOFS * l_numberOfBuffers +
                                                             NUMBER_OF_ALIGNED_DERS * l_numberOfDerivatives;

      // jump over region
      l_offset += m_meshStructure[l_cluster].numberOfCopyRegionCells[l_region];
    }
  }
}
#endif

void seissol::initializers::MemoryManager::allocateConstantData() {
  // allocate cell local data
  LocalIntegrationData       *l_local       = (LocalIntegrationData*)       m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells) * sizeof( LocalIntegrationData ),       1, MEMKIND_CONSTANT );
  NeighboringIntegrationData *l_neighboring = (NeighboringIntegrationData*) m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells) * sizeof( NeighboringIntegrationData ), 1, MEMKIND_CONSTANT );
  CellMaterialData           *material      = (CellMaterialData*)           m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells) * sizeof( CellMaterialData ),           1, MEMKIND_CONSTANT );
#ifdef USE_PLASTICITY
  PlasticityData             *plasticity    = (PlasticityData*)             m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells) * sizeof( PlasticityData ),             1, MEMKIND_CONSTANT );
#endif // USE_PLASTICITY

  // store per-cluster locations of the data
#ifdef USE_MPI
  m_copyCellData     = (struct CellData*) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( struct CellData ), 1, MEMKIND_CONSTANT );
#endif
  m_interiorCellData = (struct CellData*) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( struct CellData ), 1, MEMKIND_CONSTANT );

  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
#ifdef USE_MPI
    // set copy cell data
    m_copyCellData[l_cluster].localIntegration       = l_local;
    m_copyCellData[l_cluster].neighboringIntegration = l_neighboring;
    m_copyCellData[l_cluster].material = material;
#ifdef USE_PLASTICITY
    m_copyCellData[l_cluster].plasticity = plasticity;
#else // USE_PLASTICITY
    m_copyCellData[l_cluster].plasticity = 0L;
#endif // USE_PLASTICITY
#endif
    // jump over copy cells
    l_local       += m_meshStructure[l_cluster].numberOfCopyCells;
    l_neighboring += m_meshStructure[l_cluster].numberOfCopyCells;
    material      += m_meshStructure[l_cluster].numberOfCopyCells;
#ifdef USE_PLASTICITY
    plasticity    += m_meshStructure[l_cluster].numberOfCopyCells;
#endif // USE_PLASTICITY

    // set interior cell data
    m_interiorCellData[l_cluster].localIntegration       = l_local;
    m_interiorCellData[l_cluster].neighboringIntegration = l_neighboring;
    m_interiorCellData[l_cluster].material               = material;
#ifdef USE_PLASTICITY
    m_interiorCellData[l_cluster].plasticity             = plasticity;
#else // USE_PLASTICITY
    m_interiorCellData[l_cluster].plasticity = 0L;
#endif // USE_PLASTICITY
    
    // jump over interior
    l_local       += m_meshStructure[l_cluster].numberOfInteriorCells;
    l_neighboring += m_meshStructure[l_cluster].numberOfInteriorCells;
    material      += m_meshStructure[l_cluster].numberOfInteriorCells;
#ifdef USE_PLASTICITY
    plasticity    += m_meshStructure[l_cluster].numberOfInteriorCells;
#endif // USE_PLASTICITY
  }
}

void seissol::initializers::MemoryManager::touchConstantData( unsigned int                i_numberOfCells,
                                                              LocalIntegrationData*       o_local,
                                                              NeighboringIntegrationData* o_neighboring ) {
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
    // zero star matrices
    for( unsigned int l_starMatrix = 0; l_starMatrix < 3; l_starMatrix++ ) {
      for( unsigned int l_entry = 0; l_entry < seissol::model::AstarT::reals; l_entry++ ) {
        o_local[l_cell].starMatrices[l_starMatrix][l_entry] = (real) 0;
      }
    }

    // zero flux solvers
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      for( unsigned int l_entry = 0; l_entry < seissol::model::AplusT::reals; l_entry++ ) {
        o_local[      l_cell].nApNm1[l_face][l_entry] = (real) 0;
      }
      for( unsigned int l_entry = 0; l_entry < seissol::model::AminusT::reals; l_entry++ ) {
        o_neighboring[l_cell].nAmNm1[l_face][l_entry] = (real) 0;
      }
    }
  }
}

void seissol::initializers::MemoryManager::initializeConstantData() {
  // iterate over clusters
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
	  // TODO touch plasticity initial loading
#ifdef USE_MPI
    touchConstantData( m_meshStructure[l_cluster].numberOfCopyCells,
                       m_copyCellData[l_cluster].localIntegration,
                       m_copyCellData[l_cluster].neighboringIntegration );
#endif
    touchConstantData( m_meshStructure[l_cluster].numberOfInteriorCells,
                       m_interiorCellData[l_cluster].localIntegration,
                       m_interiorCellData[l_cluster].neighboringIntegration );

  }
}


void seissol::initializers::MemoryManager::allocateInternalState() {
  // derive sizes
#ifdef USE_MPI
  unsigned int l_ghostSize    = 0;
  unsigned int l_copySize     = 0;
#endif // USE_MPI
  unsigned int l_interiorSize = 0;
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {

#ifdef USE_MPI
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      l_ghostSize    += NUMBER_OF_ALIGNED_DOFS * m_numberOfGhostRegionBuffers[l_cluster][l_region];
      l_ghostSize    += NUMBER_OF_ALIGNED_DERS * m_numberOfGhostRegionDerivatives[l_cluster][l_region];

      l_copySize     += NUMBER_OF_ALIGNED_DOFS * m_numberOfCopyRegionBuffers[l_cluster][l_region];
      l_copySize     += NUMBER_OF_ALIGNED_DERS * m_numberOfCopyRegionDerivatives[l_cluster][l_region];
    }
#endif // USE_MPI

    l_interiorSize += NUMBER_OF_ALIGNED_DOFS * m_numberOfInteriorBuffers[l_cluster];
    l_interiorSize += NUMBER_OF_ALIGNED_DERS * m_numberOfInteriorDerivatives[l_cluster];
  }

#ifdef USE_MPI
  m_internalState.ghostLayer   = (real*) m_memoryAllocator.allocateMemory( l_ghostSize    * sizeof( real ),
                                                                           PAGESIZE_HEAP,
                                                                           MEMKIND_TIMEDOFS                );

  m_internalState.copyLayer    = (real*) m_memoryAllocator.allocateMemory( l_copySize     * sizeof( real ),
                                                                           PAGESIZE_HEAP,
                                                                           MEMKIND_TIMEDOFS                );
#endif // USE_MPI

  m_internalState.interiorTime = (real*) m_memoryAllocator.allocateMemory( l_interiorSize * sizeof( real ),
                                                                           PAGESIZE_HEAP,
                                                                           MEMKIND_TIMEDOFS                );

  /*
   * buffers / derivatives / face neighbors
   */
  m_internalState.buffers       = (real**)      m_memoryAllocator.allocateMemory( m_totalNumberOfCells * sizeof( real*    ), 1, MEMKIND_TIMEDOFS );
  m_internalState.derivatives   = (real**)      m_memoryAllocator.allocateMemory( m_totalNumberOfCells * sizeof( real*    ), 1, MEMKIND_TIMEDOFS );
  m_internalState.faceNeighbors = (real*(*)[4]) m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells) * sizeof( real*[4] ), 1, MEMKIND_TIMEDOFS );

  /*
   * dofs
   */
  m_internalState.dofs         = (real(*)[NUMBER_OF_ALIGNED_DOFS]) m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells)*sizeof( real[NUMBER_OF_ALIGNED_DOFS] ),
                                                                                                     PAGESIZE_HEAP,
                                                                                                     MEMKIND_DOFS
                                                                                                   );

#ifdef USE_PLASTICITY
  // TODO pstrain is updated every time step but only used for output -> maybe useful candidate for high bandwidth memory
  m_internalState.pstrain      = (real(*)[7]) m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells)*sizeof( real[7] ),
		  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 PAGESIZE_HEAP,
																				 seissol::memory::Standard);

  m_internalState.Energy = (real(*)[3]) m_memoryAllocator.allocateMemory( (m_totalNumberOfCopyCells + m_totalNumberOfInteriorCells)*sizeof(real[3]),
	  	  	  	  	  	  	  	                                               PAGESIZE_HEAP,
			                                                                   seissol::memory::Standard);

#else // USE_PLASTICITY
  m_internalState.pstrain = 0L;
  m_internalState.Energy = 0L;
#endif // USE_PLASTICITY
}

void seissol::initializers::MemoryManager::initializeFaceNeighbors() {
  unsigned int l_cell = 0;
  unsigned int l_ghostOffset = 0;

  // iterate over clusters
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    l_ghostOffset += m_meshStructure[l_cluster].numberOfGhostCells;

    // iterate over copy layer and interior
    for( unsigned int l_clusterCell = 0; l_clusterCell < m_meshStructure[l_cluster].numberOfCopyCells + m_meshStructure[l_cluster].numberOfInteriorCells; l_clusterCell++ ) {
      // get cell information for the current cell
      CellLocalInformation *l_cellInformation = 0L; // TODO check if NULL pointer is correct in serial version

      if( l_clusterCell < m_meshStructure[l_cluster].numberOfCopyCells ) {
#ifdef USE_MPI
        l_cellInformation = m_copyCellInformation[l_cluster] + l_clusterCell;
#endif // USE_MPI
      }
      else {
        l_cellInformation = m_interiorCellInformation[l_cluster] + l_clusterCell - m_meshStructure[l_cluster].numberOfCopyCells;
      }

      // set up pointers to face neighbors.
      for( int l_face = 0; l_face < 4; l_face++ ) {
        // default case: face neighbor is a regular element
        if( l_cellInformation->faceTypes[l_face] == regular
         || l_cellInformation->faceTypes[l_face] == periodic
         || l_cellInformation->faceTypes[l_face] == dynamicRupture ) {
          // neighboring cell provides derivatives
          if( (l_cellInformation->ltsSetup >> l_face) % 2 ) {
            m_internalState.faceNeighbors[l_cell][l_face] = m_internalState.derivatives[ l_cellInformation->faceNeighborIds[l_face] ];
          }
          // neighboring cell provides a time buffer
          else {
            m_internalState.faceNeighbors[l_cell][l_face] = m_internalState.buffers[ l_cellInformation->faceNeighborIds[l_face] ];
          }
          assert( m_internalState.faceNeighbors[l_cell][l_face] != NULL );
        }
        // free surface boundary
        else if( l_cellInformation->faceTypes[l_face] == freeSurface ) {
          if( (l_cellInformation->ltsSetup >> l_face) % 2 == 0 ) { // free surface on buffers
            m_internalState.faceNeighbors[l_cell][l_face] = m_internalState.buffers[l_cell+l_ghostOffset];
          }
          else { // free surface on derivatives
            m_internalState.faceNeighbors[l_cell][l_face] = m_internalState.derivatives[l_cell+l_ghostOffset];
          }
          assert( m_internalState.faceNeighbors[l_cell][l_face] != NULL );
        }
        // absorbing
        else if( l_cellInformation->faceTypes[l_face] == outflow ) {
          // NULL pointer; absorbing: data is not used
          m_internalState.faceNeighbors[l_cell][l_face] = NULL;
        }
        else {
          // assert all cases are covered
          assert( false );
        }
      }
      l_cell++;
    }
  }
}

void seissol::initializers::MemoryManager::initializeInternalState() {
  // offset in cells
  unsigned int l_offset = 0;

  // pointers to the clusters
#ifdef USE_MPI
  real *l_ghostPointer    = m_internalState.ghostLayer;
  real *l_copyPointer     = m_internalState.copyLayer;
#endif
  real *l_interiorPointer = m_internalState.interiorTime;

  // initialize the pointers of the internal state
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
#ifdef USE_MPI
    /*
     * ghost layer
     */
    InternalState::setUpLayerPointers( ghost,
                                       m_meshStructure[l_cluster].numberOfRegions,
                                       m_meshStructure[l_cluster].numberOfGhostRegionCells,
                                       m_ghostCellInformation[l_cluster],
                                       m_numberOfGhostRegionBuffers[l_cluster],
                                       m_numberOfGhostRegionDerivatives[l_cluster],
                                       l_ghostPointer,
                                       m_internalState.buffers     + l_offset,
                                       m_internalState.derivatives + l_offset );

    // update the offset
    l_offset  += m_meshStructure[l_cluster].numberOfGhostCells;

    // update the ghost pointer
    l_ghostPointer += NUMBER_OF_ALIGNED_DOFS * m_numberOfGhostBuffers[     l_cluster ];
    l_ghostPointer += NUMBER_OF_ALIGNED_DERS * m_numberOfGhostDerivatives[ l_cluster ];

    /*
     * Copy layer
     */
    InternalState::setUpLayerPointers( copy,
                                       m_meshStructure[l_cluster].numberOfRegions,
                                       m_meshStructure[l_cluster].numberOfCopyRegionCells,
                                       m_copyCellInformation[l_cluster],
                                       m_numberOfCopyRegionBuffers[l_cluster],
                                       m_numberOfCopyRegionDerivatives[l_cluster],
                                       l_copyPointer,
                                       m_internalState.buffers     + l_offset,
                                       m_internalState.derivatives + l_offset );

    // update the offset
    l_offset += m_meshStructure[l_cluster].numberOfCopyCells;

    // update the copy pointer
    l_copyPointer += NUMBER_OF_ALIGNED_DOFS * m_numberOfCopyBuffers[    l_cluster];
    l_copyPointer += NUMBER_OF_ALIGNED_DERS * m_numberOfCopyDerivatives[l_cluster];
#endif

    /*
     * Interior
     */
    InternalState::setUpInteriorPointers( m_meshStructure[l_cluster].numberOfInteriorCells,
                                          m_interiorCellInformation[l_cluster],
                                          m_numberOfInteriorBuffers[l_cluster],
                                          m_numberOfInteriorDerivatives[l_cluster],
                                          l_interiorPointer,
                                          m_internalState.buffers     + l_offset,
                                          m_internalState.derivatives + l_offset  );

    // update the offset
    l_offset          += m_meshStructure[l_cluster].numberOfInteriorCells;

    // update the cluster pointer
    l_interiorPointer += NUMBER_OF_ALIGNED_DOFS * m_numberOfInteriorBuffers[l_cluster];
    l_interiorPointer += NUMBER_OF_ALIGNED_DERS * m_numberOfInteriorDerivatives[l_cluster];
  }

  // initialize face neighbors
  initializeFaceNeighbors();
}

void seissol::initializers::MemoryManager::allocateCells() {
  // allocate cells struct per time cluster
  m_cells = (Cells*) m_memoryAllocator.allocateMemory( m_numberOfClusters * sizeof( Cells ), 1 );
}

void seissol::initializers::MemoryManager::touchDofs( unsigned int   i_numberOfCells,
                                                      real         (*o_dofs)[NUMBER_OF_ALIGNED_DOFS]  ) {
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
    for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
      // zero dofs
      o_dofs[l_cell][l_dof]  = (real) 0;
    }
  }
}

void seissol::initializers::MemoryManager::touchTime( unsigned int   i_numberOfCells,
                                                      real         **o_buffers,
                                                      real         **o_derivatives ) {
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
    // touch buffers
    if( o_buffers[l_cell] != NULL ) {
      for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
          // zero time integration buffers
          o_buffers[l_cell][l_dof] = (real) 0;
      }
    }

    // touch derivatives
    if( o_derivatives[l_cell] != NULL ) {
      for( unsigned int l_derivative = 0; l_derivative < NUMBER_OF_ALIGNED_DERS; l_derivative++ ) {
        o_derivatives[l_cell][l_derivative] = (real) 0;
      }
    }
  }
}

void seissol::initializers::MemoryManager::touchPstrain(unsigned int   i_numberOfCells,
                                                      	real         (*o_pstrain)[7] )
{
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
    for( unsigned int l_dof = 0; l_dof < 7; l_dof++ ) {
      // zero pstrain output
      o_pstrain[l_cell][l_dof]  = (real) 0;
    }
  }
}

void seissol::initializers::MemoryManager::touchEnergy(unsigned int   i_numberOfCells,
                                                       real  (*o_Energy)[3] )
{
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
	  for( unsigned int l_var = 0; l_var < 2; l_var++ ) {
      // zero plastic energy output
      o_Energy[l_cell][l_var] = (real) 0;
  }
 }
}

void seissol::initializers::MemoryManager::initializeCells() {
  /*
   * Pointers to dofs
   */
  real (*l_dofsPointer)[NUMBER_OF_ALIGNED_DOFS] = m_internalState.dofs;
  real (*pstrainPointer)[7] = m_internalState.pstrain;
  real (*EnergyPointer)[3] = m_internalState.Energy;

  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
#ifdef USE_MPI
    m_cells[l_cluster].copyDofs      = l_dofsPointer;
#ifdef USE_PLASTICITY
    m_cells[l_cluster].copyPstrain   = pstrainPointer;
    m_cells[l_cluster].copyEnergy   = EnergyPointer;
#else // USE_PLASTICITY
    m_cells[l_cluster].copyPstrain   = 0L;
    m_cells[l_cluster].copyEnergy   = 0L;
#endif // USE_PLASTICITY
#endif // USE_MPI

    // jump over copy layer
    l_dofsPointer += m_meshStructure[l_cluster].numberOfCopyCells;
    pstrainPointer += m_meshStructure[l_cluster].numberOfCopyCells;
    EnergyPointer += m_meshStructure[l_cluster].numberOfCopyCells;

    // set interior dofs
    m_cells[l_cluster].interiorDofs = l_dofsPointer;
#ifdef USE_PLASTICITY
    m_cells[l_cluster].interiorPstrain   = pstrainPointer;
    m_cells[l_cluster].interiorEnergy   = EnergyPointer;
#else // USE_PLASTICITY
    m_cells[l_cluster].interiorPstrain   = 0L;
    m_cells[l_cluster].interiorEnergy   = 0L;
#endif // USE_PLASTICITY

    // jump over interior
    l_dofsPointer += m_meshStructure[l_cluster].numberOfInteriorCells;
    pstrainPointer += m_meshStructure[l_cluster].numberOfInteriorCells;
    EnergyPointer += m_meshStructure[l_cluster].numberOfInteriorCells;
  }

  /*
   * Time and face pointers
   */
  // offset in cells
  unsigned int l_totalOffset = 0;
  unsigned int l_copyInteriorOffset = 0;

  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
    // jump over ghost layer
    l_totalOffset += m_meshStructure[l_cluster].numberOfGhostCells;

#ifdef USE_MPI
    // set copy layer
    m_cells[l_cluster].copyBuffers           = m_internalState.buffers       + l_totalOffset;
    m_cells[l_cluster].copyDerivatives       = m_internalState.derivatives   + l_totalOffset;
    m_cells[l_cluster].copyFaceNeighbors     = m_internalState.faceNeighbors + l_copyInteriorOffset;
#endif // USE_MPI

    // jump over copy layer
    l_totalOffset        += m_meshStructure[l_cluster].numberOfCopyCells;
    l_copyInteriorOffset += m_meshStructure[l_cluster].numberOfCopyCells;

    // set interior
    m_cells[l_cluster].interiorBuffers       = m_internalState.buffers       + l_totalOffset;
    m_cells[l_cluster].interiorDerivatives   = m_internalState.derivatives   + l_totalOffset;
    m_cells[l_cluster].interiorFaceNeighbors = m_internalState.faceNeighbors + l_copyInteriorOffset;

    // jump of over interior
    l_totalOffset        += m_meshStructure[l_cluster].numberOfInteriorCells;
    l_copyInteriorOffset += m_meshStructure[l_cluster].numberOfInteriorCells;
  }

  /*
   * Touch memory
   */
  for( unsigned int l_cluster = 0; l_cluster < m_numberOfClusters; l_cluster++ ) {
#ifdef USE_MPI
      // initialize copy dofs
      touchDofs( m_meshStructure[l_cluster].numberOfCopyCells,
                 m_cells[l_cluster].copyDofs );

      // initialize copy buffers / derivatives
      touchTime( m_meshStructure[l_cluster].numberOfCopyCells,
                 m_cells[l_cluster].copyBuffers,
                 m_cells[l_cluster].copyDerivatives );

#ifdef USE_PLASTICITY
      touchPstrain( m_meshStructure[l_cluster].numberOfCopyCells,
    		  m_cells[l_cluster].copyPstrain );
      touchEnergy( m_meshStructure[l_cluster].numberOfCopyCells,
    		  m_cells[l_cluster].copyEnergy );
#endif // USE_PLASTICITY
#endif

      // initialize interior dofs
      touchDofs( m_meshStructure[l_cluster].numberOfInteriorCells,
                 m_cells[l_cluster].interiorDofs );

      // initialize interior buffers / derivatives
      touchTime( m_meshStructure[l_cluster].numberOfInteriorCells,
                 m_cells[l_cluster].interiorBuffers,
                 m_cells[l_cluster].interiorDerivatives );
#ifdef USE_PLASTICITY
      touchPstrain( m_meshStructure[l_cluster].numberOfInteriorCells,
    		  m_cells[l_cluster].interiorPstrain );
      touchEnergy( m_meshStructure[l_cluster].numberOfInteriorCells,
    		  m_cells[l_cluster].interiorEnergy );
#endif // USE_PLASTICITY

  }
}

void seissol::initializers::MemoryManager::initializeMemoryLayout( struct TimeStepping         &i_timeStepping,
                                                                   struct MeshStructure        *i_meshStructure,
                                                                   struct CellLocalInformation *io_cellLocalInformation ) {
  // store mesh structure and the number of time clusters
  m_meshStructure = i_meshStructure;
  m_numberOfClusters = i_timeStepping.numberOfLocalClusters;

  // correct LTS-information in the ghost layer
  correctGhostRegionSetups( io_cellLocalInformation );

  // set up the layers
  setUpLayers( io_cellLocalInformation );

  // derive the layouts of the layers
  deriveLayerLayouts();

  // allocate memory for the constant data
  allocateConstantData();

  // initialize the contant data
  initializeConstantData();

  // allocate the internal state
  allocateInternalState();

  // initialize the internal state
  initializeInternalState();

  // allocate the cells
  allocateCells();

  // initialize the cells
  initializeCells();

#ifdef USE_MPI
  // initialize the communication structure
  initializeCommunicationStructure();
#endif
}

void seissol::initializers::MemoryManager::getMemoryLayout( unsigned int                    i_cluster,
                                                            struct MeshStructure          *&o_meshStructure,
#ifdef USE_MPI
                                                            struct CellLocalInformation   *&o_copyCellInformation,
#endif
                                                            struct CellLocalInformation   *&o_interiorCellInformation,
                                                            struct GlobalData             *&o_globalData,
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
                                                            struct GlobalData             *&o_globalDataCopies,
#endif
#ifdef USE_MPI
                                                            struct CellData               *&o_copyCellData,
#endif
                                                            struct CellData               *&o_interiorCellData,
                                                            struct Cells                  *&o_cells ) {
  o_meshStructure           =  m_meshStructure + i_cluster;
#ifdef USE_MPI
  o_copyCellInformation     =  m_copyCellInformation[i_cluster];
#endif
  o_interiorCellInformation =  m_interiorCellInformation[i_cluster];
  o_globalData              = &m_globalData;
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
  o_globalDataCopies        =  m_globalDataCopies;
#endif
#ifdef USE_MPI
  o_copyCellData            =  m_copyCellData     + i_cluster;
#endif
  o_interiorCellData        =  m_interiorCellData + i_cluster;
  o_cells                   = &m_cells[i_cluster];
}
