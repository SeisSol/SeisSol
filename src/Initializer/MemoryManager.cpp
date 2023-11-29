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
#include "GlobalData.h"

#include <Kernels/common.hpp>
#include <Kernels/Touch.h>
#include <unordered_set>
#include <cmath>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include "device.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#endif // ACL_DEVICE


void seissol::initializers::MemoryManager::initialize()
{
  // initialize global matrices
  GlobalDataInitializerOnHost::init(m_globalDataOnHost, m_memoryAllocator, MEMKIND_GLOBAL);
  if constexpr (seissol::isDeviceOn()) {
    // the serial order for initialization is needed for some (older) driver versions on some GPUs
    bool serialize = false;
    const char* envvalue = std::getenv("SEISSOL_SERIAL_NODE_DEVICE_INIT");
    if (envvalue != nullptr) {
      if (strcmp(envvalue, "1") == 0) {
        serialize = true;
      }
      else if (strcmp(envvalue, "0") == 0) {
        serialize = false;
      }
      else {
        logError() << "Invalid value for \"SEISSOL_SERIAL_NODE_DEVICE_INIT\"";
      }
    }
    if (serialize) {
      logInfo(MPI::mpi.rank()) << "Initializing device global data on a node in serial order.";
      MPI::mpi.serialOrderExecute([&]() {
        GlobalDataInitializerOnDevice::init(m_globalDataOnDevice, m_memoryAllocator, memory::DeviceGlobalMemory);
      }, MPI::mpi.sharedMemComm());
    }
    else {
      GlobalDataInitializerOnDevice::init(m_globalDataOnDevice, m_memoryAllocator, memory::DeviceGlobalMemory);
    }
  }
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
      m_meshStructure[tc].ghostRegionSizes[l_region] = tensor::Q::size() * l_numberOfBuffers +
                                                       yateto::computeFamilySize<tensor::dQ>() * l_numberOfDerivatives;

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
      m_meshStructure[tc].copyRegionSizes[l_region] = tensor::Q::size() * l_numberOfBuffers +
                                                      yateto::computeFamilySize<tensor::dQ>() * l_numberOfDerivatives;

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
      if (cellInformation[cell].faceTypes[face] == FaceType::regular ||
	  cellInformation[cell].faceTypes[face] == FaceType::periodic ||
	  cellInformation[cell].faceTypes[face] == FaceType::dynamicRupture) {
        // neighboring cell provides derivatives
        if( (cellInformation[cell].ltsSetup >> face) % 2 ) {
          faceNeighbors[cell][face] = derivatives[ cellInformation[cell].faceNeighborIds[face] ];
        }
        // neighboring cell provides a time buffer
        else {
          faceNeighbors[cell][face] = buffers[ cellInformation[cell].faceNeighborIds[face] ];
        }
        assert(faceNeighbors[cell][face] != nullptr);
      }
      // boundaries using local cells
      else if (cellInformation[cell].faceTypes[face] == FaceType::freeSurface ||
	       cellInformation[cell].faceTypes[face] == FaceType::freeSurfaceGravity ||
	       cellInformation[cell].faceTypes[face] == FaceType::dirichlet ||
	       cellInformation[cell].faceTypes[face] == FaceType::analytical) {
        if( (cellInformation[cell].ltsSetup >> face) % 2 == 0 ) { // free surface on buffers
          faceNeighbors[cell][face] = layer.var(m_lts.buffers)[cell];
        }
        else { // free surface on derivatives
          faceNeighbors[cell][face] = layer.var(m_lts.derivatives)[cell];
        }
        assert(faceNeighbors[cell][face] != nullptr);
      }
      // absorbing
      else if( cellInformation[cell].faceTypes[face] == FaceType::outflow ) {
        // NULL pointer; absorbing: data is not used
        faceNeighbors[cell][face] = nullptr;
      }
      else {
        // assert all cases are covered
        assert(false);
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

void seissol::initializers::MemoryManager::fixateLtsTree(struct TimeStepping& i_timeStepping,
                                                         struct MeshStructure*i_meshStructure,
                                                         unsigned* numberOfDRCopyFaces,
                                                         unsigned* numberOfDRInteriorFaces,
                                                         bool usePlasticity) {
  // store mesh structure and the number of time clusters
  m_meshStructure = i_meshStructure;

  // Setup tree variables
  m_lts.addTo(m_ltsTree, usePlasticity);
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
  m_dynRup->addTo(m_dynRupTree);

  m_dynRupTree.setNumberOfTimeClusters(i_timeStepping.numberOfGlobalClusters);
  m_dynRupTree.fixate();

  for (unsigned tc = 0; tc < m_dynRupTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_dynRupTree.child(tc);
    cluster.child<Ghost>().setNumberOfCells(0);
    if (tc >= i_timeStepping.numberOfLocalClusters) {
        cluster.child<Copy>().setNumberOfCells(0);
        cluster.child<Interior>().setNumberOfCells(0);
    } else {
        cluster.child<Copy>().setNumberOfCells(numberOfDRCopyFaces[tc]);
        cluster.child<Interior>().setNumberOfCells(numberOfDRInteriorFaces[tc]);
    }
  }

  m_dynRupTree.allocateVariables();
  m_dynRupTree.touchVariables();

#ifdef ACL_DEVICE
  MemoryManager::deriveRequiredScratchpadMemoryForDr(m_dynRupTree, *m_dynRup.get());
  m_dynRupTree.allocateScratchPads();
#endif
}

void seissol::initializers::MemoryManager::fixateBoundaryLtsTree() {
  seissol::initializers::LayerMask ghostMask(Ghost);

  // Boundary face tree
  m_boundary.addTo(m_boundaryTree);
  m_boundaryTree.setNumberOfTimeClusters(m_ltsTree.numChildren());
  m_boundaryTree.fixate();

  // First count the number of faces with relevant boundary condition.
  for (unsigned tc = 0; tc < m_boundaryTree.numChildren(); ++tc) {
    auto& cluster = m_boundaryTree.child(tc);
    cluster.child<Ghost>().setNumberOfCells(0);
    cluster.child<Copy>().setNumberOfCells(0);
    cluster.child<Interior>().setNumberOfCells(0);
  }

  // Iterate over layers of standard lts tree and face lts tree together.
  auto layer = m_ltsTree.beginLeaf(ghostMask), boundaryLayer = m_boundaryTree.beginLeaf(ghostMask);
    for (;
       layer != m_ltsTree.endLeaf() && boundaryLayer != m_boundaryTree.endLeaf();
       ++layer, ++boundaryLayer) {
    CellLocalInformation* cellInformation = layer->var(m_lts.cellInformation);

    unsigned numberOfBoundaryFaces = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          ++numberOfBoundaryFaces;
        }
      }
    }
    boundaryLayer->setNumberOfCells(numberOfBoundaryFaces);
  }
  m_boundaryTree.allocateVariables();
  m_boundaryTree.touchVariables();

  // The boundary tree is now allocated, now we only need to map from cell lts
  // to face lts.
  // We do this by, once again, iterating over both trees at the same time.
  for (auto layer = m_ltsTree.beginLeaf(ghostMask), boundaryLayer = m_boundaryTree.beginLeaf(ghostMask);
       layer != m_ltsTree.endLeaf() && boundaryLayer != m_boundaryTree.endLeaf();
       ++layer, ++boundaryLayer) {
    auto* cellInformation = layer->var(m_lts.cellInformation);
    auto* boundaryMapping = layer->var(m_lts.boundaryMapping);
    auto* faceInformation = boundaryLayer->var(m_boundary.faceInformation);

    auto boundaryFace = 0;
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          boundaryMapping[cell][face].nodes = faceInformation[boundaryFace].nodes;
          boundaryMapping[cell][face].TData = faceInformation[boundaryFace].TData;
          boundaryMapping[cell][face].TinvData = faceInformation[boundaryFace].TinvData;
          boundaryMapping[cell][face].easiBoundaryMap = faceInformation[boundaryFace].easiBoundaryMap;
          boundaryMapping[cell][face].easiBoundaryConstant = faceInformation[boundaryFace].easiBoundaryConstant;
          ++boundaryFace;
        } else {
          boundaryMapping[cell][face].nodes = nullptr;
          boundaryMapping[cell][face].TData = nullptr;
          boundaryMapping[cell][face].TinvData = nullptr;
          boundaryMapping[cell][face].easiBoundaryMap = nullptr;
          boundaryMapping[cell][face].easiBoundaryConstant = nullptr;
        }
      }
    }
  }
}

void seissol::initializers::MemoryManager::deriveFaceDisplacementsBucket()
{
  for (auto layer = m_ltsTree.beginLeaf(m_lts.faceDisplacements.mask); layer != m_ltsTree.endLeaf(); ++layer) {
    CellLocalInformation* cellInformation = layer->var(m_lts.cellInformation);
    real* (*displacements)[4] = layer->var(m_lts.faceDisplacements);
    CellMaterialData* cellMaterialData = layer->var(m_lts.material);

    unsigned numberOfFaces = 0;
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      for (unsigned int face = 0; face < 4; ++face) {
        if (requiresDisplacement(cellInformation[cell],
                                 cellMaterialData[cell],
                                 face)) {
          // We add the base address later when the bucket is allocated
          // +1 is necessary as we want to reserve the nullptr for cell without displacement.
          // Thanks to this hack, the array contains a constant plus the offset of the current
          // cell.
          displacements[cell][face] =
              static_cast<real*>(nullptr) + 1 + numberOfFaces * tensor::faceDisplacement::size();
          ++numberOfFaces;
        } else {
          displacements[cell][face] = nullptr;
        }
      }
    }
    layer->setBucketSize(m_lts.faceDisplacementsBuffer, numberOfFaces * 1 * tensor::faceDisplacement::size() * sizeof(real));
  }
}

#ifdef ACL_DEVICE
void seissol::initializers::MemoryManager::deriveRequiredScratchpadMemoryForWp(LTSTree& ltsTree, LTS& lts) {
  constexpr size_t totalDerivativesSize = yateto::computeFamilySize<tensor::dQ>();
  constexpr size_t nodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto layer = ltsTree.beginLeaf(Ghost); layer != ltsTree.endLeaf(); ++layer) {

    CellLocalInformation *cellInformation = layer->var(lts.cellInformation);
    std::unordered_set<real *> registry{};
    real *(*faceNeighbors)[4] = layer->var(lts.faceNeighbors);

    unsigned derivativesCounter{0};
    unsigned integratedDofsCounter{0};
    unsigned nodalDisplacementsCounter{0};

    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      bool needsScratchMemForDerivatives = (cellInformation[cell].ltsSetup >> 9) % 2 == 0;
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++integratedDofsCounter;

      // include data provided by ghost layers
      for (unsigned face = 0; face < 4; ++face) {
        real *neighbourBuffer = faceNeighbors[cell][face];

        // check whether a neighbour element idofs has not been counted twice
        if ((registry.find(neighbourBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighbourBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] != FaceType::outflow &&
                cellInformation[cell].faceTypes[face] != FaceType::dynamicRupture) {

              bool isNeighbProvidesDerivatives = ((cellInformation[cell].ltsSetup >> face) % 2) == 1;
              if (isNeighbProvidesDerivatives) {
                ++integratedDofsCounter;
              }
              registry.insert(neighbourBuffer);
            }
          }
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::freeSurfaceGravity) {
          ++nodalDisplacementsCounter;
        }

      }
    }
    layer->setScratchpadSize(lts.integratedDofsScratch,
                             integratedDofsCounter * tensor::I::size() * sizeof(real));
    layer->setScratchpadSize(lts.derivativesScratch,
                             derivativesCounter * totalDerivativesSize * sizeof(real));
    layer->setScratchpadSize(lts.nodalAvgDisplacements,
                             nodalDisplacementsCounter * nodalDisplacementsSize * sizeof(real));
  }
}

void seissol::initializers::MemoryManager::deriveRequiredScratchpadMemoryForDr(
    LTSTree &ltsTree,
    DynamicRupture& dynRup) {
  constexpr size_t idofsSize = tensor::Q::size() * sizeof(real);
  for (auto layer = ltsTree.beginLeaf(); layer != ltsTree.endLeaf(); ++layer) {
    const auto layerSize = layer->getNumberOfCells();
    layer->setScratchpadSize(dynRup.idofsPlusOnDevice, idofsSize * layerSize);
    layer->setScratchpadSize(dynRup.idofsMinusOnDevice, idofsSize * layerSize);
  }
}
#endif

void seissol::initializers::MemoryManager::initializeFaceDisplacements()
{
  for (auto layer = m_ltsTree.beginLeaf(m_lts.faceDisplacements.mask); layer != m_ltsTree.endLeaf(); ++layer) {
    if (layer->getBucketSize(m_lts.faceDisplacementsBuffer) == 0) {
      continue;
    }
    real* (*displacements)[4] = layer->var(m_lts.faceDisplacements);
    real* bucket = static_cast<real*>(layer->bucket(m_lts.faceDisplacementsBuffer));

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(layer, displacements, bucket)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (displacements[cell][face] != nullptr) {
          // Remove constant part that was added in deriveDisplacementsBucket.
          // We then have the pointer offset that needs to be added to the bucket.
          // The final value of this pointer then points to a valid memory address
          // somewhere in the bucket.
          displacements[cell][face] = bucket + ((displacements[cell][face] - static_cast<real*>(nullptr)) - 1);
          for (unsigned dof = 0; dof < tensor::faceDisplacement::size(); ++dof) {
            // zero displacements
            displacements[cell][face][dof] = static_cast<real>(0.0);
          }
        }
      }
    }
  }
}

void seissol::initializers::MemoryManager::initializeMemoryLayout()
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
      l_ghostSize    += sizeof(real) * tensor::Q::size() * m_numberOfGhostRegionBuffers[tc][l_region];
      l_ghostSize    += sizeof(real) * yateto::computeFamilySize<tensor::dQ>() * m_numberOfGhostRegionDerivatives[tc][l_region];

      l_copySize     += sizeof(real) * tensor::Q::size() * m_numberOfCopyRegionBuffers[tc][l_region];
      l_copySize     += sizeof(real) * yateto::computeFamilySize<tensor::dQ>() * m_numberOfCopyRegionDerivatives[tc][l_region];
    }
#endif // USE_MPI
    l_interiorSize += sizeof(real) * tensor::Q::size() * m_numberOfInteriorBuffers[tc];
    l_interiorSize += sizeof(real) * yateto::computeFamilySize<tensor::dQ>() * m_numberOfInteriorDerivatives[tc];

    cluster.child<Ghost>().setBucketSize(m_lts.buffersDerivatives, l_ghostSize);
    cluster.child<Copy>().setBucketSize(m_lts.buffersDerivatives, l_copySize);
    cluster.child<Interior>().setBucketSize(m_lts.buffersDerivatives, l_interiorSize);
  }

  deriveFaceDisplacementsBucket();

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

#ifdef ACL_DEVICE
  void* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
  for (auto it = m_ltsTree.beginLeaf(); it != m_ltsTree.endLeaf(); ++it) {
    if (it->getBucketSize(m_lts.buffersDerivatives) > 0) {
      void* data = it->bucket(m_lts.buffersDerivatives);
      device::DeviceInstance::getInstance().algorithms.touchMemory(
        reinterpret_cast<real*>(data),
        it->getBucketSize(m_lts.buffersDerivatives) / sizeof(real),
        true, stream);
    }
  }
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#else
  for (auto it = m_ltsTree.beginLeaf(); it != m_ltsTree.endLeaf(); ++it) {
    real** buffers = it->var(m_lts.buffers);
    real** derivatives = it->var(m_lts.derivatives);
    kernels::touchBuffersDerivatives(buffers, derivatives, it->getNumberOfCells());
  }
#endif

#ifdef USE_MPI
  // initialize the communication structure
  initializeCommunicationStructure();
#endif

  initializeFaceDisplacements();

#ifdef ACL_DEVICE
  seissol::initializers::MemoryManager::deriveRequiredScratchpadMemoryForWp(m_ltsTree, m_lts);
  m_ltsTree.allocateScratchPads();
#endif
}

std::pair<MeshStructure *, CompoundGlobalData>
seissol::initializers::MemoryManager::getMemoryLayout(unsigned int i_cluster) {
  MeshStructure *meshStructure = m_meshStructure + i_cluster;

  CompoundGlobalData globalData{};
  globalData.onHost = &m_globalDataOnHost;
  globalData.onDevice = nullptr;
  if constexpr (seissol::isDeviceOn()) {
    globalData.onDevice = &m_globalDataOnDevice;
  }

  return std::make_pair(meshStructure, globalData);
}

void seissol::initializers::MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (fileNameStr != "") {
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}


#ifdef ACL_DEVICE
void seissol::initializers::MemoryManager::recordExecutionPaths(bool usePlasticity) {
  recording::CompositeRecorder<seissol::initializers::LTS> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

  if (usePlasticity) {
    recorder.addRecorder(new recording::PlasticityRecorder);
  }

  for (LTSTree::leaf_iterator it = m_ltsTree.beginLeaf(Ghost); it != m_ltsTree.endLeaf(); ++it) {
    recorder.record(m_lts, *it);
  }

  recording::CompositeRecorder<seissol::initializers::DynamicRupture> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (LTSTree::leaf_iterator it = m_dynRupTree.beginLeaf(Ghost); it != m_dynRupTree.endLeaf(); ++it) {
    drRecorder.record(*m_dynRup, *it);
  }
}
#endif // ACL_DEVICE

bool seissol::initializers::isAcousticSideOfElasticAcousticInterface(CellMaterialData &material,
                                              unsigned int face) {
#ifdef USE_ANISOTROPIC
  return false;
#else
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face].mu > eps && material.local.mu < eps;
#endif
}
bool seissol::initializers::isElasticSideOfElasticAcousticInterface(CellMaterialData &material,
                                             unsigned int face) {
#ifdef USE_ANISOTROPIC
  return false;
#else
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.local.mu > eps && material.neighbor[face].mu < eps;
#endif
}

bool seissol::initializers::isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
#ifndef USE_ANISOTROPIC
  return isAcousticSideOfElasticAcousticInterface(material, face) || isElasticSideOfElasticAcousticInterface(material, face);
#else
  return false;
#endif
}


bool seissol::initializers::requiresDisplacement(CellLocalInformation cellLocalInformation,
                                                 CellMaterialData &material,
                                                 unsigned int face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::freeSurface
  || faceType == FaceType::freeSurfaceGravity
  || isAtElasticAcousticInterface(material, face);
}

bool seissol::initializers::requiresNodalFlux(FaceType f) {
  return (f == FaceType::freeSurfaceGravity
          || f == FaceType::dirichlet
          || f == FaceType::analytical);
}

void seissol::initializers::MemoryManager::initializeFrictionLaw() {
  const int rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Initialize Friction Model";

  const auto factory = seissol::dr::factory::getFactory(m_dynRupParameters);
  auto product = factory->produce();
  m_dynRup = std::move(product.ltsTree);
  m_DRInitializer = std::move(product.initializer);
  m_FrictionLaw = std::move(product.frictionLaw);
  m_faultOutputManager = std::move(product.output);
}

void seissol::initializers::MemoryManager::initFaultOutputManager(const std::string& backupTimeStamp) {
  // TODO: switch m_dynRup to shared or weak pointer
  if (m_dynRupParameters->isDynamicRuptureEnabled) {
    m_faultOutputManager->setInputParam(*m_inputParams, seissol::SeisSol::main.meshReader());
    m_faultOutputManager->setLtsData(&m_ltsTree,
                                     &m_lts,
                                     &m_ltsLut,
                                     &m_dynRupTree,
                                     m_dynRup.get());
    m_faultOutputManager->setBackupTimeStamp(backupTimeStamp);
    m_faultOutputManager->init();

  }
}


void seissol::initializers::MemoryManager::initFrictionData() {
  if (m_dynRupParameters->isDynamicRuptureEnabled) {

    m_DRInitializer->initializeFault(m_dynRup.get(), &m_dynRupTree);

#ifdef ACL_DEVICE
    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(m_FrictionLaw.get())) {
      impl->initSyclQueue();

      LayerMask mask = seissol::initializers::LayerMask(Ghost);
      auto maxSize = m_dynRupTree.getMaxClusterSize(mask);
      impl->setMaxClusterSize(maxSize);

      impl->allocateAuxiliaryMemory();
      impl->copyStaticDataToDevice();
    }
#endif // ACL_DEVICE
  }
}

