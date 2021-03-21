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
#include <yateto.h>

#include <Kernels/common.hpp>
#include <generated_code/tensor.h>
#include <unordered_set>
#include <cmath>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include <Solver/Pipeline/DrPipeline.h>
#endif //ACL_DEVICE

void seissol::initializers::MemoryManager::initialize()
{
  // initialize global matrices
  GlobalDataInitializerOnHost::init(m_globalDataOnHost, m_memoryAllocator, MEMKIND_GLOBAL);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(m_globalDataOnDevice, m_memoryAllocator, memory::DeviceGlobalMemory);
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
      for (unsigned dof = 0; dof < tensor::Q::size(); ++dof) {
          // zero time integration buffers
          buffer[dof] = (real) 0;
      }
    }

    // touch derivatives
    real* derivative = derivatives[cell];
    if (derivative != NULL) {
      for (unsigned dof = 0; dof < yateto::computeFamilySize<tensor::dQ>(); ++dof ) {
        derivative[dof] = (real) 0;
      }
    }
  }
}

void seissol::initializers::MemoryManager::fixateLtsTree(struct TimeStepping& i_timeStepping,
                                                         struct MeshStructure*i_meshStructure,
                                                         unsigned* numberOfDRCopyFaces,
                                                         unsigned* numberOfDRInteriorFaces) {
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

  for (unsigned tc = 0; tc < m_dynRupTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_dynRupTree.child(tc);
    cluster.child<Ghost>().setNumberOfCells(0);
    cluster.child<Copy>().setNumberOfCells(numberOfDRCopyFaces[tc]);
    cluster.child<Interior>().setNumberOfCells(numberOfDRInteriorFaces[tc]);
  }

  m_dynRupTree.allocateVariables();
  m_dynRupTree.touchVariables();

#ifdef ACL_DEVICE
  constexpr size_t QInterpolatedSize = CONVERGENCE_ORDER * tensor::QInterpolated::size() * sizeof(real);
  constexpr size_t imposedStateSize = tensor::QInterpolated::size() * sizeof(real);
  constexpr size_t idofsSize = tensor::Q::size() * sizeof(real);
  for (auto layer = m_dynRupTree.beginLeaf(); layer != m_dynRupTree.endLeaf(); ++layer) {
    const auto layerSize = layer->getNumberOfCells();
    layer->setScratchpadSize(m_dynRup.QInterpolatedPlusOnDevice, QInterpolatedSize * layerSize);
    layer->setScratchpadSize(m_dynRup.QInterpolatedMinusOnDevice, QInterpolatedSize * layerSize);
    layer->setScratchpadSize(m_dynRup.idofsPlusOnDevice, idofsSize * layerSize);
    layer->setScratchpadSize(m_dynRup.idofsMinusOnDevice, idofsSize * layerSize);

    constexpr auto UpperStageFactor = dr::pipeline::DrPipeline::TailSize * dr::pipeline::DrPipeline::DefaultBatchSize;
    constexpr auto LowerStageFactor = dr::pipeline::DrPipeline::NumStages * dr::pipeline::DrPipeline::DefaultBatchSize;
    layer->setScratchpadSize(m_dynRup.QInterpolatedPlusOnHost, UpperStageFactor * QInterpolatedSize);
    layer->setScratchpadSize(m_dynRup.QInterpolatedMinusOnHost, UpperStageFactor * QInterpolatedSize);
    layer->setScratchpadSize(m_dynRup.imposedStatePlusOnHost, LowerStageFactor * imposedStateSize);
    layer->setScratchpadSize(m_dynRup.imposedStateMinusOnHost, LowerStageFactor *  imposedStateSize);
  }
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

void seissol::initializers::MemoryManager::deriveDisplacementsBucket()
{
  for (auto layer = m_ltsTree.beginLeaf(m_lts.displacements.mask); layer != m_ltsTree.endLeaf(); ++layer) {
    CellLocalInformation* cellInformation = layer->var(m_lts.cellInformation);
    real** displacements = layer->var(m_lts.displacements);
    CellMaterialData* cellMaterialData = layer->var(m_lts.material);

    unsigned numberOfCells = 0;
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      bool hasFreeSurface = false;
      for (unsigned int face = 0; face < 4; ++face) {
        const auto faceType = cellInformation[cell].faceTypes[face];
        hasFreeSurface = hasFreeSurface
                         || faceType == FaceType::freeSurface
                         || faceType == FaceType::freeSurfaceGravity
                         || isAtElasticAcousticInterface(cellMaterialData[cell], face);
      }
      if (hasFreeSurface) {
        // We add the base address later when the bucket is allocated
        // +1 is necessary as we want to reserve the nullptr for cell without displacement.
        // Thanks to this hack, the array contains a constant plus the offset of the current
        // cell.
        displacements[cell] = static_cast<real*>(nullptr) + 1 + numberOfCells * tensor::displacement::size();
        ++numberOfCells;
      } else {
        displacements[cell] = nullptr;
      }
    }
    layer->setBucketSize(m_lts.displacementsBuffer, numberOfCells * tensor::displacement::size() * sizeof(real));
  }
}

#ifdef ACL_DEVICE
void seissol::initializers::MemoryManager::deriveRequiredScratchpadMemory() {
  constexpr size_t totalDerivativesSize = yateto::computeFamilySize<tensor::dQ>();
  kernels::NeighborData::Loader loader;

  for (auto layer = m_ltsTree.beginLeaf(Ghost); layer != m_ltsTree.endLeaf(); ++layer) {

    CellLocalInformation *cellInformation = layer->var(m_lts.cellInformation);
    std::unordered_set<real *> registry{};
    real *(*faceNeighbors)[4] = layer->var(m_lts.faceNeighbors);

    unsigned derivativesCounter{0};
    unsigned idofsCounter{0};

    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {

      bool needsScratchMemForDerivatives = (cellInformation[cell].ltsSetup >> 9) % 2 == 0;
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++idofsCounter;

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
                ++idofsCounter;
              }
              registry.insert(neighbourBuffer);
            }
          }
        }
      }
    }
    layer->setScratchpadSize(m_lts.idofsScratch,
                             idofsCounter * tensor::I::size() * sizeof(real));
    layer->setScratchpadSize(m_lts.derivativesScratch,
                             derivativesCounter * totalDerivativesSize * sizeof(real));
  }
}
#endif

void seissol::initializers::MemoryManager::initializeDisplacements()
{
  for (auto layer = m_ltsTree.beginLeaf(m_lts.displacements.mask); layer != m_ltsTree.endLeaf(); ++layer) {
    if (layer->getBucketSize(m_lts.displacementsBuffer) == 0) {
      continue;
    }
    real** displacements = layer->var(m_lts.displacements);
    real* bucket = static_cast<real*>(layer->bucket(m_lts.displacementsBuffer));

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
      if (displacements[cell] != nullptr) {
        // Remove constant part that was added in deriveDisplacementsBucket.
        // We then have the pointer offset that needs to be added to the bucket.
        // The final value of this pointer then points to a valid memory address
        // somewhere in the bucket.
        displacements[cell] = bucket + ((displacements[cell] - static_cast<real*>(nullptr)) - 1);
        for (unsigned dof = 0; dof < tensor::displacement::size(); ++dof) {
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

  deriveDisplacementsBucket();

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

  for (auto it = m_ltsTree.beginLeaf(); it != m_ltsTree.endLeaf(); ++it) {
    touchBuffersDerivatives(*it);
  }

#ifdef USE_MPI
  // initialize the communication structure
  initializeCommunicationStructure();
#endif

  initializeDisplacements();

#ifdef ACL_DEVICE
  deriveRequiredScratchpadMemory();
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
    std::cout << "initializeEasiBoundaryReader with file: " << fileName << std::endl;
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}


#ifdef ACL_DEVICE
void seissol::initializers::MemoryManager::recordExecutionPaths() {
  recording::CompositeRecorder<seissol::initializers::LTS> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

#ifdef USE_PLASTICITY
  recorder.addRecorder(new recording::PlasticityRecorder);
#endif

  for (LTSTree::leaf_iterator it = m_ltsTree.beginLeaf(Ghost); it != m_ltsTree.endLeaf(); ++it) {
    recorder.record(m_lts, *it);
  }

  recording::CompositeRecorder<seissol::initializers::DynamicRupture> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (LTSTree::leaf_iterator it = m_dynRupTree.beginLeaf(Ghost); it != m_dynRupTree.endLeaf(); ++it) {
    drRecorder.record(m_dynRup, *it);
  }
}
#endif // ACL_DEVICE


bool seissol::initializers::isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
  //TODO (LK): implement more general coupling
#ifndef USE_ANISOTROPIC
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.local.mu > eps && material.neighbor[face].mu < eps;
#else
  return false;
#endif
}

bool seissol::initializers::requiresNodalFlux(FaceType f) {
  return (f == FaceType::freeSurfaceGravity
          || f == FaceType::dirichlet
          || f == FaceType::analytical);
}
