// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Memory/MemoryAllocator.h"
#include "SeisSol.h"
#include "MemoryManager.h"
#include "InternalState.h"
#include "Memory/Tree/Layer.h"
#include <cstddef>
#include <yateto.h>
#include <unordered_set>
#include <cmath>
#include <type_traits>
#include "Memory/GlobalData.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/Common.h"
#include "Kernels/Touch.h"

#include <DynamicRupture/Misc.h>

#include "Common/Iterator.h"

#include "generated_code/tensor.h"

#ifdef ACL_DEVICE
#include "BatchRecorders/Recorders.h"
#include "device.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#endif // ACL_DEVICE

void seissol::initializer::MemoryManager::initialize()
{
  // initialize global matrices
  GlobalDataInitializerOnHost::init(m_globalDataOnHost, m_memoryAllocator, memory::Standard);
  if constexpr (seissol::isDeviceOn()) {
    GlobalDataInitializerOnDevice::init(m_globalDataOnDevice, m_memoryAllocator, memory::DeviceGlobalMemory);
  }
}

void seissol::initializer::MemoryManager::correctGhostRegionSetups()
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

void seissol::initializer::MemoryManager::deriveLayerLayouts() {
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
void seissol::initializer::MemoryManager::initializeCommunicationStructure() {
  // reset mpi requests
  for( unsigned int l_cluster = 0; l_cluster < m_ltsTree.numChildren(); l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      m_meshStructure[l_cluster].sendRequests[l_region] = MPI_REQUEST_NULL;
      m_meshStructure[l_cluster].receiveRequests[l_region] = MPI_REQUEST_NULL;
    }
  }

#ifdef ACL_DEVICE
  const auto allocationPlace = seissol::initializer::AllocationPlace::Device;
#else
  const auto allocationPlace = seissol::initializer::AllocationPlace::Host;
#endif

  /*
   * ghost layer
   */
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
    real* ghostStart = static_cast<real*>(cluster.child<Ghost>().bucket(m_lts.buffersDerivatives, allocationPlace));
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
#ifdef ACL_DEVICE
    real** buffers = copy.var(m_lts.buffersDevice);
    real** derivatives = copy.var(m_lts.derivativesDevice);
#else
    real** buffers = copy.var(m_lts.buffers);
    real** derivatives = copy.var(m_lts.derivatives);
#endif
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

void seissol::initializer::MemoryManager::initializeFaceNeighbors( unsigned    cluster,
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
#ifdef ACL_DEVICE
  real** buffersDevice = m_ltsTree.var(m_lts.buffersDevice);          // faceNeighborIds are ltsIds and not layer-local
  real** derivativesDevice = m_ltsTree.var(m_lts.derivativesDevice);  // faceNeighborIds are ltsIds and not layer-local
  real *(*faceNeighborsDevice)[4] = layer.var(m_lts.faceNeighborsDevice);
#endif
  auto* cellInformation = layer.var(m_lts.cellInformation);
  auto* secondaryInformation = layer.var(m_lts.secondaryInformation);

  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    for (unsigned face = 0; face < 4; ++face) {
      if (cellInformation[cell].faceTypes[face] == FaceType::Regular ||
	  cellInformation[cell].faceTypes[face] == FaceType::Periodic ||
	  cellInformation[cell].faceTypes[face] == FaceType::DynamicRupture) {
        // neighboring cell provides derivatives
        if( (cellInformation[cell].ltsSetup >> face) % 2 ) {
          faceNeighbors[cell][face] = derivatives[ secondaryInformation[cell].faceNeighborIds[face] ];
#ifdef ACL_DEVICE
          faceNeighborsDevice[cell][face] = derivativesDevice[ secondaryInformation[cell].faceNeighborIds[face] ];
#endif
        }
        // neighboring cell provides a time buffer
        else {
          faceNeighbors[cell][face] = buffers[ secondaryInformation[cell].faceNeighborIds[face] ];
#ifdef ACL_DEVICE
          faceNeighborsDevice[cell][face] = buffersDevice[ secondaryInformation[cell].faceNeighborIds[face] ];
#endif
        }
        assert(faceNeighbors[cell][face] != nullptr);
      }
      // boundaries using local cells
      else if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurface ||
	       cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity ||
	       cellInformation[cell].faceTypes[face] == FaceType::Dirichlet ||
	       cellInformation[cell].faceTypes[face] == FaceType::Analytical) {
        if( (cellInformation[cell].ltsSetup >> face) % 2 == 0 ) { // free surface on buffers
          faceNeighbors[cell][face] = layer.var(m_lts.buffers)[cell];
#ifdef ACL_DEVICE
          faceNeighborsDevice[cell][face] = layer.var(m_lts.buffersDevice)[cell];
#endif
        }
        else { // free surface on derivatives
          faceNeighbors[cell][face] = layer.var(m_lts.derivatives)[cell];
#ifdef ACL_DEVICE
          faceNeighborsDevice[cell][face] = layer.var(m_lts.derivativesDevice)[cell];
#endif
        }
        assert(faceNeighbors[cell][face] != nullptr);
      }
      // absorbing
      else if( cellInformation[cell].faceTypes[face] == FaceType::Outflow ) {
        // NULL pointer; absorbing: data is not used
        faceNeighbors[cell][face] = nullptr;
#ifdef ACL_DEVICE
        faceNeighborsDevice[cell][face] = nullptr;
#endif
      }
      else {
        // assert all cases are covered
        assert(false);
      }
    }
  }
}

void seissol::initializer::MemoryManager::initializeBuffersDerivatives() {
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

#ifdef ACL_DEVICE
    #ifdef USE_MPI
    /*
     * ghost layer
     */
    InternalState::setUpLayerPointers( m_meshStructure[tc].numberOfRegions,
                                       m_meshStructure[tc].numberOfGhostRegionCells,
                                       cluster.child<Ghost>().var(m_lts.cellInformation),
                                       m_numberOfGhostRegionBuffers[tc],
                                       m_numberOfGhostRegionDerivatives[tc],
                                       static_cast<real*>(cluster.child<Ghost>().bucket(m_lts.buffersDerivatives, seissol::initializer::AllocationPlace::Device)),
                                       cluster.child<Ghost>().var(m_lts.buffersDevice),
                                       cluster.child<Ghost>().var(m_lts.derivativesDevice) );

    /*
     * Copy layer
     */
    InternalState::setUpLayerPointers( m_meshStructure[tc].numberOfRegions,
                                       m_meshStructure[tc].numberOfCopyRegionCells,
                                       cluster.child<Copy>().var(m_lts.cellInformation),
                                       m_numberOfCopyRegionBuffers[tc],
                                       m_numberOfCopyRegionDerivatives[tc],
                                       static_cast<real*>(cluster.child<Copy>().bucket(m_lts.buffersDerivatives, seissol::initializer::AllocationPlace::Device)),
                                       cluster.child<Copy>().var(m_lts.buffersDevice),
                                       cluster.child<Copy>().var(m_lts.derivativesDevice) );
#endif

    /*
     * Interior
     */
    InternalState::setUpInteriorPointers( m_meshStructure[tc].numberOfInteriorCells,
                                          cluster.child<Interior>().var(m_lts.cellInformation),
                                          m_numberOfInteriorBuffers[tc],
                                          m_numberOfInteriorDerivatives[tc],
                                          static_cast<real*>(cluster.child<Interior>().bucket(m_lts.buffersDerivatives, seissol::initializer::AllocationPlace::Device)),
                                          cluster.child<Interior>().var(m_lts.buffersDevice),
                                          cluster.child<Interior>().var(m_lts.derivativesDevice)  );
#endif
  }
}

void seissol::initializer::MemoryManager::fixateLtsTree(struct TimeStepping& i_timeStepping,
                                                         struct MeshStructure*i_meshStructure,
                                                         unsigned* numberOfDRCopyFaces,
                                                         unsigned* numberOfDRInteriorFaces,
                                                         bool usePlasticity) {
  // store mesh structure and the number of time clusters
  m_meshStructure = i_meshStructure;

  // Setup tree variables
  m_lts.addTo(m_ltsTree, usePlasticity);
  seissolInstance.postProcessor().allocateMemory(&m_ltsTree);
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

void seissol::initializer::MemoryManager::fixateBoundaryLtsTree() {
  seissol::initializer::LayerMask ghostMask(Ghost);

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
  for (auto [layer, boundaryLayer] : seissol::common::zip(m_ltsTree.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);

    unsigned numberOfBoundaryFaces = 0;
    const auto layerSize = layer.getNumberOfCells();
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layerSize; ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          ++numberOfBoundaryFaces;
        }
      }
    }
    boundaryLayer.setNumberOfCells(numberOfBoundaryFaces);
  }
  m_boundaryTree.allocateVariables();
  m_boundaryTree.touchVariables();

  // The boundary tree is now allocated, now we only need to map from cell lts
  // to face lts.
  // We do this by, once again, iterating over both trees at the same time.
  for (auto [layer, boundaryLayer] : seissol::common::zip(m_ltsTree.leaves(ghostMask), m_boundaryTree.leaves(ghostMask))) {
    auto* cellInformation = layer.var(m_lts.cellInformation);
    auto* boundaryMapping = layer.var(m_lts.boundaryMapping);
    auto* boundaryMappingDevice = layer.var(m_lts.boundaryMappingDevice);
    auto* faceInformation = boundaryLayer.var(m_boundary.faceInformation, AllocationPlace::Host);
    auto* faceInformationDevice = boundaryLayer.var(m_boundary.faceInformation, AllocationPlace::Device);

    auto boundaryFace = 0;
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          boundaryMapping[cell][face].nodes = faceInformation[boundaryFace].nodes;
          boundaryMapping[cell][face].TData = faceInformation[boundaryFace].TData;
          boundaryMapping[cell][face].TinvData = faceInformation[boundaryFace].TinvData;
          boundaryMapping[cell][face].easiBoundaryMap = faceInformation[boundaryFace].easiBoundaryMap;
          boundaryMapping[cell][face].easiBoundaryConstant = faceInformation[boundaryFace].easiBoundaryConstant;
          boundaryMappingDevice[cell][face].nodes = faceInformationDevice[boundaryFace].nodes;
          boundaryMappingDevice[cell][face].TData = faceInformationDevice[boundaryFace].TData;
          boundaryMappingDevice[cell][face].TinvData = faceInformationDevice[boundaryFace].TinvData;
          boundaryMappingDevice[cell][face].easiBoundaryMap = faceInformationDevice[boundaryFace].easiBoundaryMap;
          boundaryMappingDevice[cell][face].easiBoundaryConstant = faceInformationDevice[boundaryFace].easiBoundaryConstant;
          ++boundaryFace;
        } else {
          boundaryMapping[cell][face].nodes = nullptr;
          boundaryMapping[cell][face].TData = nullptr;
          boundaryMapping[cell][face].TinvData = nullptr;
          boundaryMapping[cell][face].easiBoundaryMap = nullptr;
          boundaryMapping[cell][face].easiBoundaryConstant = nullptr;
          boundaryMappingDevice[cell][face].nodes = nullptr;
          boundaryMappingDevice[cell][face].TData = nullptr;
          boundaryMappingDevice[cell][face].TinvData = nullptr;
          boundaryMappingDevice[cell][face].easiBoundaryMap = nullptr;
          boundaryMappingDevice[cell][face].easiBoundaryConstant = nullptr;
        }
      }
    }
  }
}

void seissol::initializer::MemoryManager::deriveFaceDisplacementsBucket()
{
  for (auto& layer : m_ltsTree.leaves(m_lts.faceDisplacements.mask)) {
    CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);
    real* (*displacements)[4] = layer.var(m_lts.faceDisplacements);
    CellMaterialData* cellMaterialData = layer.var(m_lts.material);

    unsigned numberOfFaces = 0;
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      for (unsigned int face = 0; face < 4; ++face) {
        if (requiresDisplacement(cellInformation[cell],
                                 cellMaterialData[cell],
                                 face)) {
          // We add the base address later when the bucket is allocated
          // +1 is necessary as we want to reserve the nullptr for cell without displacement.
          // Thanks to this hack, the array contains a constant plus the offset of the current
          // cell.
          displacements[cell][face] =
              reinterpret_cast<real*>(1 + numberOfFaces * tensor::faceDisplacement::size());
          ++numberOfFaces;
        } else {
          displacements[cell][face] = nullptr;
        }
      }
    }
    layer.setBucketSize(m_lts.faceDisplacementsBuffer, numberOfFaces * 1 * tensor::faceDisplacement::size() * sizeof(real));
  }
}

#ifdef ACL_DEVICE
void seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(LTSTree& ltsTree, LTS& lts) {
  constexpr size_t totalDerivativesSize = yateto::computeFamilySize<tensor::dQ>();
  constexpr size_t nodalDisplacementsSize = tensor::averageNormalDisplacement::size();

  for (auto& layer : ltsTree.leaves(Ghost)) {

    CellLocalInformation *cellInformation = layer.var(lts.cellInformation);
    std::unordered_set<real *> registry{};
    real *(*faceNeighbors)[4] = layer.var(lts.faceNeighborsDevice);

    std::size_t derivativesCounter{0};
    std::size_t integratedDofsCounter{0};
    std::size_t nodalDisplacementsCounter{0};
    std::size_t analyticCounter = 0;

    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      bool needsScratchMemForDerivatives = (cellInformation[cell].ltsSetup >> 9) % 2 == 0;
      if (needsScratchMemForDerivatives) {
        ++derivativesCounter;
      }
      ++integratedDofsCounter;

      // include data provided by ghost layers
      for (unsigned face = 0; face < 4; ++face) {
        real *neighborBuffer = faceNeighbors[cell][face];

        // check whether a neighbour element idofs has not been counted twice
        if ((registry.find(neighborBuffer) == registry.end())) {

          // maybe, because of BCs, a pointer can be a nullptr, i.e. skip it
          if (neighborBuffer != nullptr) {
            if (cellInformation[cell].faceTypes[face] != FaceType::Outflow &&
                cellInformation[cell].faceTypes[face] != FaceType::DynamicRupture) {

              bool isNeighbProvidesDerivatives = ((cellInformation[cell].ltsSetup >> face) % 2) == 1;
              if (isNeighbProvidesDerivatives) {
                ++integratedDofsCounter;
              }
              registry.insert(neighborBuffer);
            }
          }
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::FreeSurfaceGravity) {
          ++nodalDisplacementsCounter;
        }

        if (cellInformation[cell].faceTypes[face] == FaceType::Analytical) {
          ++analyticCounter;
        }
      }
    }
    layer.setScratchpadSize(lts.integratedDofsScratch,
                             integratedDofsCounter * tensor::I::size() * sizeof(real));
    layer.setScratchpadSize(lts.derivativesScratch,
                             derivativesCounter * totalDerivativesSize * sizeof(real));
    layer.setScratchpadSize(lts.nodalAvgDisplacements,
                             nodalDisplacementsCounter * nodalDisplacementsSize * sizeof(real));
    layer.setScratchpadSize(lts.analyticScratch,
                             analyticCounter * tensor::INodal::size() * sizeof(real));
  }
}

void seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForDr(
    LTSTree &ltsTree,
    DynamicRupture& dynRup) {
  constexpr size_t idofsSize = tensor::Q::size() * sizeof(real);
  for (auto& layer : ltsTree.leaves()) {
    const auto layerSize = layer.getNumberOfCells();
    layer.setScratchpadSize(dynRup.idofsPlusOnDevice, idofsSize * layerSize);
    layer.setScratchpadSize(dynRup.idofsMinusOnDevice, idofsSize * layerSize);
  }
}
#endif

void seissol::initializer::MemoryManager::initializeFaceDisplacements()
{
  for (auto& layer : m_ltsTree.leaves(m_lts.faceDisplacements.mask)) {
    if (layer.getBucketSize(m_lts.faceDisplacementsBuffer) == 0) {
      continue;
    }
    real* (*displacements)[4] = layer.var(m_lts.faceDisplacements);
    real* bucket = static_cast<real*>(layer.bucket(m_lts.faceDisplacementsBuffer));
    real* (*displacementsDevice)[4] = layer.var(m_lts.faceDisplacementsDevice);
    real* bucketDevice = static_cast<real*>(layer.bucket(m_lts.faceDisplacementsBuffer, seissol::initializer::AllocationPlace::Device));

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(layer, displacements, bucket, displacementsDevice, bucketDevice)
#endif // _OPENMP
    for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        if (displacements[cell][face] != nullptr) {
          // Remove constant part that was added in deriveDisplacementsBucket.
          // We then have the pointer offset that needs to be added to the bucket.
          // The final value of this pointer then points to a valid memory address
          // somewhere in the bucket.
          auto offset = (reinterpret_cast<std::intptr_t>(displacements[cell][face]) - 1);
          displacements[cell][face] = bucket + offset;
          displacementsDevice[cell][face] = bucketDevice + offset;
          for (unsigned dof = 0; dof < tensor::faceDisplacement::size(); ++dof) {
            // zero displacements
            displacements[cell][face][dof] = static_cast<real>(0.0);
          }
        }
      }
    }
  }
}

void seissol::initializer::MemoryManager::initializeMemoryLayout()
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
  for (auto& layer : m_ltsTree.leaves()) {
    if (layer.getBucketSize(m_lts.buffersDerivatives) > 0) {
      void* data = layer.bucket(m_lts.buffersDerivatives, seissol::initializer::AllocationPlace::Device);
      device::DeviceInstance::getInstance().algorithms.touchMemory(
        reinterpret_cast<real*>(data),
        layer.getBucketSize(m_lts.buffersDerivatives) / sizeof(real),
        true, stream);
    }
  }
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
  for (auto& layer : m_ltsTree.leaves()) {
    real** buffers = layer.var(m_lts.buffers);
    real** derivatives = layer.var(m_lts.derivatives);
    kernels::touchBuffersDerivatives(buffers, derivatives, layer.getNumberOfCells());
  }

#ifdef USE_MPI
  // initialize the communication structure
  initializeCommunicationStructure();
#endif

  initializeFaceDisplacements();

#ifdef ACL_DEVICE
  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(m_ltsTree, m_lts);
  m_ltsTree.allocateScratchPads();
#endif
}

std::pair<MeshStructure *, CompoundGlobalData>
seissol::initializer::MemoryManager::getMemoryLayout(unsigned int i_cluster) {
  MeshStructure *meshStructure = m_meshStructure + i_cluster;

  CompoundGlobalData globalData{};
  globalData.onHost = &m_globalDataOnHost;
  globalData.onDevice = nullptr;
  if constexpr (seissol::isDeviceOn()) {
    globalData.onDevice = &m_globalDataOnDevice;
  }

  return std::make_pair(meshStructure, globalData);
}

void seissol::initializer::MemoryManager::initializeEasiBoundaryReader(const char* fileName) {
  const auto fileNameStr = std::string{fileName};
  if (fileNameStr != "") {
    m_easiBoundary = EasiBoundary(fileNameStr);
  }
}


#ifdef ACL_DEVICE
void seissol::initializer::MemoryManager::recordExecutionPaths(bool usePlasticity) {
  recording::CompositeRecorder<seissol::initializer::LTS> recorder;
  recorder.addRecorder(new recording::LocalIntegrationRecorder);
  recorder.addRecorder(new recording::NeighIntegrationRecorder);

  if (usePlasticity) {
    recorder.addRecorder(new recording::PlasticityRecorder);
  }

  for (auto& layer : m_ltsTree.leaves(Ghost)) {
    recorder.record(m_lts, layer);
  }

  recording::CompositeRecorder<seissol::initializer::DynamicRupture> drRecorder;
  drRecorder.addRecorder(new recording::DynamicRuptureRecorder);
  for (auto& layer : m_dynRupTree.leaves(Ghost)) {
    drRecorder.record(*m_dynRup, layer);
  }
}
#endif // ACL_DEVICE

bool seissol::initializer::isAcousticSideOfElasticAcousticInterface(CellMaterialData &material,
                                              unsigned int face) {
#ifdef USE_ANISOTROPIC
  return false;
#else
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face].getMuBar() > eps && material.local.getMuBar() < eps;
#endif
}
bool seissol::initializer::isElasticSideOfElasticAcousticInterface(CellMaterialData &material,
                                             unsigned int face) {
#ifdef USE_ANISOTROPIC
  return false;
#else
  constexpr auto eps = std::numeric_limits<real>::epsilon();
  return material.local.getMuBar() > eps && material.neighbor[face].getMuBar() < eps;
#endif
}

bool seissol::initializer::isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
#ifndef USE_ANISOTROPIC
  return isAcousticSideOfElasticAcousticInterface(material, face) || isElasticSideOfElasticAcousticInterface(material, face);
#else
  return false;
#endif
}


bool seissol::initializer::requiresDisplacement(CellLocalInformation cellLocalInformation,
                                                 CellMaterialData &material,
                                                 unsigned int face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::FreeSurface
  || faceType == FaceType::FreeSurfaceGravity
  || isAtElasticAcousticInterface(material, face);
}

bool seissol::initializer::requiresNodalFlux(FaceType f) {
  return (f == FaceType::FreeSurfaceGravity
          || f == FaceType::Dirichlet
          || f == FaceType::Analytical);
}

void seissol::initializer::MemoryManager::initializeFrictionLaw() {
  const auto drParameters = std::make_shared<seissol::initializer::parameters::DRParameters>(m_seissolParams->drParameters);
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(drParameters->frictionLawType).c_str()
    << "(" << static_cast<int>(drParameters->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (drParameters->isThermalPressureOn ? "on" : "off");

  const auto factory = seissol::dr::factory::getFactory(drParameters, seissolInstance);
  auto product = factory->produce();
  m_dynRup = std::move(product.ltsTree);
  m_DRInitializer = std::move(product.initializer);
  m_FrictionLaw = std::move(product.frictionLaw);
  m_FrictionLawDevice = std::move(product.frictionLawDevice);
  m_faultOutputManager = std::move(product.output);
}

void seissol::initializer::MemoryManager::initFaultOutputManager(const std::string& backupTimeStamp) {
  // TODO: switch m_dynRup to shared or weak pointer
  if (m_seissolParams->drParameters.isDynamicRuptureEnabled) {
    m_faultOutputManager->setInputParam(seissolInstance.meshReader());
    m_faultOutputManager->setLtsData(&m_ltsTree,
                                     &m_lts,
                                     &m_ltsLut,
                                     &m_dynRupTree,
                                     m_dynRup.get());
    m_faultOutputManager->setBackupTimeStamp(backupTimeStamp);
    m_faultOutputManager->init();

  }
}


void seissol::initializer::MemoryManager::initFrictionData() {
  if (m_seissolParams->drParameters.isDynamicRuptureEnabled) {

    m_DRInitializer->initializeFault(m_dynRup.get(), &m_dynRupTree);

#ifdef ACL_DEVICE
    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(m_FrictionLawDevice.get())) {

      const auto mask = seissol::initializer::LayerMask(Ghost);

      impl->allocateAuxiliaryMemory();
    }
#endif // ACL_DEVICE
  }
}

void seissol::initializer::MemoryManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  if (place == seissol::initializer::AllocationPlace::Device) {
    logInfo() << "Synchronizing data... (host->device)";
  }
  else {
    logInfo() << "Synchronizing data... (device->host)";
  }
  const auto& defaultStream = device::DeviceInstance::getInstance().api->getDefaultStream();
  m_ltsTree.synchronizeTo(place, defaultStream);
  m_dynRupTree.synchronizeTo(place, defaultStream);
  m_boundaryTree.synchronizeTo(place, defaultStream);
  device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
#endif
}

