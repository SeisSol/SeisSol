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

#include "MemoryManager.h"
#include "InternalState.h"
#include "Tree/Layer.h"
#include <cstddef>
#include <yateto.h>
#include <unordered_set>
#include <cmath>
#include <type_traits>
#include "GlobalData.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Kernels/Common.h"
#include "Kernels/Touch.h"

#include <DynamicRupture/Misc.h>

#include "Common/Iterator.h"

#include "generated_code/tensor.h"
#include "GlobalData.h"
#include "InternalState.h"
#include "SeisSol.h"


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

    unsigned int lOffset = 0;
    for( unsigned int lRegion = 0; lRegion < m_meshStructure[tc].numberOfRegions; lRegion++ ) {
      // iterate over ghost cells
      for( unsigned int lCell = 0; lCell < m_meshStructure[tc].numberOfGhostRegionCells[lRegion]; lCell++ ) {
        if( lCell < m_meshStructure[tc].numberOfGhostRegionDerivatives[lRegion] ) {
          // assert the cell provides derivatives
          assert( (cellInformation[lOffset+lCell].ltsSetup >> 9)%2 );

          // reset possible buffers
          cellInformation[lOffset+lCell].ltsSetup &= ( ~(1 << 8 ) );
          cellInformation[lOffset+lCell].ltsSetup &= ( ~(1 << 10) );
        } else {
          // assert the cell provides buffers
          assert( (cellInformation[lOffset+lCell].ltsSetup >> 8)%2 );

          // reset possible derivatives
          cellInformation[lOffset+lCell].ltsSetup &= ( ~(1 << 9 ) );
        }
      }
      // update offset with ghost region size
      lOffset +=  m_meshStructure[tc].numberOfGhostRegionCells[lRegion];
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
    unsigned int lGhostOffset = 0;
    unsigned int lCopyOffset  = 0;
    for( unsigned int lRegion = 0; lRegion < m_meshStructure[tc].numberOfRegions; lRegion++ ) {
      m_numberOfGhostRegionBuffers[     tc][lRegion] = 0;
      m_numberOfGhostRegionDerivatives[ tc][lRegion] = 0;

      m_numberOfCopyRegionBuffers[      tc][lRegion] = 0;
      m_numberOfCopyRegionDerivatives[  tc][lRegion] = 0;

      // iterate over all cells of this clusters ghost layer
      for( unsigned int lCell = 0; lCell < m_meshStructure[tc].numberOfGhostRegionCells[lRegion]; lCell++ ) {
        // ensure that either buffers or derivatives are used; not both!
        bool lBuffer      = ( ghostCellInformation[lCell+lGhostOffset].ltsSetup >> 8 ) % 2;
        bool lDerivatives = ( ghostCellInformation[lCell+lGhostOffset].ltsSetup >> 9 ) % 2;

        if( (lBuffer && lDerivatives) || ( lBuffer || lDerivatives ) == false ) logError() << "invalid ghost lts setup" << lBuffer << lDerivatives;

        // check if this cell requires a buffer and/or derivatives
        if( ( ghostCellInformation[lCell+lGhostOffset].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfGhostRegionBuffers[    tc][lRegion]++;
        if( ( ghostCellInformation[lCell+lGhostOffset].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfGhostRegionDerivatives[tc][lRegion]++;
      }
      lGhostOffset += m_meshStructure[tc].numberOfGhostRegionCells[lRegion];

      // iterate over all cells of this clusters copy layer
      for( unsigned int lCell = 0; lCell < m_meshStructure[tc].numberOfCopyRegionCells[lRegion]; lCell++ ) {
        // assert that buffers or derivatives are requested
        assert( ( ( copyCellInformation[lCell+lCopyOffset].ltsSetup >> 8 ) % 2 ||
                  ( copyCellInformation[lCell+lCopyOffset].ltsSetup >> 9 ) % 2 )
                == true );

        // check if this cell requires a buffer and/or derivatives
        if( ( copyCellInformation[lCell+lCopyOffset].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfCopyRegionBuffers[    tc][lRegion]++;
        if( ( copyCellInformation[lCell+lCopyOffset].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfCopyRegionDerivatives[tc][lRegion]++;
      }
      lCopyOffset += m_meshStructure[tc].numberOfCopyRegionCells[lRegion];

      // update number of buffers and derivatives in the ghost and copy layer
      m_numberOfGhostBuffers[     tc ]  += m_numberOfGhostRegionBuffers[     tc][lRegion];
      m_numberOfGhostDerivatives[ tc ]  += m_numberOfGhostRegionDerivatives[ tc][lRegion];
      m_numberOfCopyBuffers[      tc ]  += m_numberOfCopyRegionBuffers[      tc][lRegion];
      m_numberOfCopyDerivatives[  tc ]  += m_numberOfCopyRegionDerivatives[  tc][lRegion];
    }
#endif // USE_MPI

    // iterate over all cells of this clusters interior
    for( unsigned int lCell = 0; lCell < m_meshStructure[tc].numberOfInteriorCells; lCell++ ) {
      // check if this cell requires a buffer and/or derivatives
      if( ( interiorCellInformation[lCell].ltsSetup >> 8 ) % 2 == 1 ) m_numberOfInteriorBuffers[    tc]++;
      if( ( interiorCellInformation[lCell].ltsSetup >> 9 ) % 2 == 1 ) m_numberOfInteriorDerivatives[tc]++;
    }
  }
}

#ifdef USE_MPI
void seissol::initializer::MemoryManager::initializeCommunicationStructure() {
  // reset mpi requests
  for( unsigned int lCluster = 0; lCluster < m_ltsTree.numChildren(); lCluster++ ) {
    for( unsigned int lRegion = 0; lRegion < m_meshStructure[lCluster].numberOfRegions; lRegion++ ) {
      m_meshStructure[lCluster].sendRequests[lRegion] = MPI_REQUEST_NULL;
      m_meshStructure[lCluster].receiveRequests[lRegion] = MPI_REQUEST_NULL;
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
      unsigned int lNumberOfDerivatives = m_meshStructure[tc].numberOfGhostRegionDerivatives[l_region];
      unsigned int lNumberOfBuffers     = m_meshStructure[tc].numberOfGhostRegionCells[l_region] - lNumberOfDerivatives;

      // set size
      m_meshStructure[tc].ghostRegionSizes[l_region] = tensor::Q::size() * lNumberOfBuffers +
                                                       yateto::computeFamilySize<tensor::dQ>() * lNumberOfDerivatives;

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
    unsigned int lOffset = 0;

    for( unsigned int lRegion = 0; lRegion < m_meshStructure[tc].numberOfRegions; lRegion++ ) {
      // derive the communication size
      unsigned int lNumberOfDerivatives = m_meshStructure[tc].numberOfCommunicatedCopyRegionDerivatives[lRegion];
      unsigned int lNumberOfBuffers     = m_meshStructure[tc].numberOfCopyRegionCells[lRegion] - lNumberOfDerivatives;

      // assert number of communicated buffers fits into total number of buffers
      assert( m_numberOfCopyRegionBuffers[tc][lRegion] >= lNumberOfBuffers );

      // set pointer to copy region start
      if( lNumberOfBuffers > 0 ) {
        m_meshStructure[tc].copyRegions[lRegion] = buffers[lNumberOfDerivatives + lOffset];
      }
      else {
        m_meshStructure[tc].copyRegions[lRegion] = derivatives[lOffset];
      }

      // assert the pointer is set
      assert( m_meshStructure[tc].copyRegions[lRegion] != NULL );

      // set size
      m_meshStructure[tc].copyRegionSizes[lRegion] = tensor::Q::size() * lNumberOfBuffers +
                                                      yateto::computeFamilySize<tensor::dQ>() * lNumberOfDerivatives;

      // jump over region
      lOffset += m_meshStructure[tc].numberOfCopyRegionCells[lRegion];
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
  CellLocalInformation* cellInformation = layer.var(m_lts.cellInformation);

  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    for (unsigned face = 0; face < 4; ++face) {
      if (cellInformation[cell].faceTypes[face] == FaceType::Regular ||
	  cellInformation[cell].faceTypes[face] == FaceType::Periodic ||
	  cellInformation[cell].faceTypes[face] == FaceType::DynamicRupture) {
        // neighboring cell provides derivatives
        if( (cellInformation[cell].ltsSetup >> face) % 2 ) {
          faceNeighbors[cell][face] = derivatives[ cellInformation[cell].faceNeighborIds[face] ];
#ifdef ACL_DEVICE
          faceNeighborsDevice[cell][face] = derivativesDevice[ cellInformation[cell].faceNeighborIds[face] ];
#endif
        }
        // neighboring cell provides a time buffer
        else {
          faceNeighbors[cell][face] = buffers[ cellInformation[cell].faceNeighborIds[face] ];
#ifdef ACL_DEVICE
          faceNeighborsDevice[cell][face] = buffersDevice[ cellInformation[cell].faceNeighborIds[face] ];
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

void seissol::initializer::MemoryManager::fixateLtsTree(struct TimeStepping& iTimeStepping,
                                                         struct MeshStructure*iMeshStructure,
                                                         unsigned* numberOfDRCopyFaces,
                                                         unsigned* numberOfDRInteriorFaces,
                                                         bool usePlasticity) {
  // store mesh structure and the number of time clusters
  m_meshStructure = iMeshStructure;

  // Setup tree variables
  m_lts.addTo(m_ltsTree, usePlasticity);
  seissolInstance.postProcessor().allocateMemory(&m_ltsTree);
  m_ltsTree.setNumberOfTimeClusters(iTimeStepping.numberOfLocalClusters);

  /// From this point, the tree layout, variables, and buckets cannot be changed anymore
  m_ltsTree.fixate();

  // Set number of cells and bucket sizes in ltstree
  for (unsigned tc = 0; tc < m_ltsTree.numChildren(); ++tc) {
    TimeCluster& cluster = m_ltsTree.child(tc);
    cluster.child<Ghost>().setNumberOfCells(iMeshStructure[tc].numberOfGhostCells);
    cluster.child<Copy>().setNumberOfCells(iMeshStructure[tc].numberOfCopyCells);
    cluster.child<Interior>().setNumberOfCells(iMeshStructure[tc].numberOfInteriorCells);
  }

  m_ltsTree.allocateVariables();
  m_ltsTree.touchVariables();

  /// Dynamic rupture tree
  for (int i = 0; i < MULTIPLE_SIMULATIONS; i++) {
    m_dynRupTree[i] = new LTSTree();
    m_dynRup[i]->addTo(*m_dynRupTree[i]);
    m_dynRupTree[i]->setNumberOfTimeClusters(iTimeStepping.numberOfGlobalClusters);
    m_dynRupTree[i]->fixate();
    for (unsigned tc = 0; tc < m_dynRupTree[i]->numChildren(); ++tc) {
      TimeCluster& cluster = m_dynRupTree[i]->child(tc);
      cluster.child<Ghost>().setNumberOfCells(0);
      if (tc >= iTimeStepping.numberOfLocalClusters) {
        cluster.child<Copy>().setNumberOfCells(0);
        cluster.child<Interior>().setNumberOfCells(0);
      } else {
        cluster.child<Copy>().setNumberOfCells(numberOfDRCopyFaces[tc]);
        cluster.child<Interior>().setNumberOfCells(numberOfDRInteriorFaces[tc]);
      }
    }
    m_dynRupTree[i]->allocateVariables();
    m_dynRupTree[i]->touchVariables();
  }

  // m_dynRup->addTo(m_dynRupTree);
  // m_dynRupTree.setNumberOfTimeClusters(i_timeStepping.numberOfGlobalClusters);
  // m_dynRupTree.fixate();

  // for (unsigned tc = 0; tc < m_dynRupTree.numChildren(); ++tc) {
  //   TimeCluster& cluster = m_dynRupTree.child(tc);
  //   cluster.child<Ghost>().setNumberOfCells(0);
  //   if (tc >= i_timeStepping.numberOfLocalClusters) {
  //       cluster.child<Copy>().setNumberOfCells(0);
  //       cluster.child<Interior>().setNumberOfCells(0);
  //   } else {
  //       cluster.child<Copy>().setNumberOfCells(numberOfDRCopyFaces[tc]);
  //       cluster.child<Interior>().setNumberOfCells(numberOfDRInteriorFaces[tc]);
  //   }
  // }

  // m_dynRupTree.allocateVariables();
  // m_dynRupTree.touchVariables();

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

    size_t lGhostSize = 0;
    size_t lCopySize = 0;
    size_t lInteriorSize = 0;
#ifdef USE_MPI
    for( unsigned int lRegion = 0; lRegion < m_meshStructure[tc].numberOfRegions; lRegion++ ) {
      lGhostSize    += sizeof(real) * tensor::Q::size() * m_numberOfGhostRegionBuffers[tc][lRegion];
      lGhostSize    += sizeof(real) * yateto::computeFamilySize<tensor::dQ>() * m_numberOfGhostRegionDerivatives[tc][lRegion];

      lCopySize     += sizeof(real) * tensor::Q::size() * m_numberOfCopyRegionBuffers[tc][lRegion];
      lCopySize     += sizeof(real) * yateto::computeFamilySize<tensor::dQ>() * m_numberOfCopyRegionDerivatives[tc][lRegion];
    }
#endif // USE_MPI
    lInteriorSize += sizeof(real) * tensor::Q::size() * m_numberOfInteriorBuffers[tc];
    lInteriorSize += sizeof(real) * yateto::computeFamilySize<tensor::dQ>() * m_numberOfInteriorDerivatives[tc];

    cluster.child<Ghost>().setBucketSize(m_lts.buffersDerivatives, lGhostSize);
    cluster.child<Copy>().setBucketSize(m_lts.buffersDerivatives, lCopySize);
    cluster.child<Interior>().setBucketSize(m_lts.buffersDerivatives, lInteriorSize);
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
seissol::initializer::MemoryManager::getMemoryLayout(unsigned int iCluster) {
  MeshStructure *meshStructure = m_meshStructure + iCluster;

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
  logInfo() << "Initialize Friction Model";

  logInfo() << "Friction law:" << dr::misc::frictionLawName(m_seissolParams->drParameters[0]->frictionLawType).c_str()
    << "(" << static_cast<int>(m_seissolParams->drParameters[0]->frictionLawType) << ")";
  logInfo() << "Thermal pressurization:" << (m_seissolParams->drParameters[0]->isThermalPressureOn ? "on" : "off");


  for (int i = 0; i < MULTIPLE_SIMULATIONS; i++) {
    auto drParameters = m_seissolParams->drParameters[i];
    const auto factory = seissol::dr::factory::getFactory(drParameters, seissolInstance, i);
    auto product = factory->produce();
    m_dynRup[i] = std::move(product.ltsTree);
    m_DRInitializer[i] = std::move(product.initializer);
    m_FrictionLaw[i] = std::move(product.frictionLaw);
    m_faultOutputManager[i] = std::move(product.output);
  }
}

void seissol::initializer::MemoryManager::initFaultOutputManager(const std::string& backupTimeStamp) {
  for (int i = 0; i < MULTIPLE_SIMULATIONS; i++) {
    if (m_seissolParams->drParameters[i]->isDynamicRuptureEnabled) {
      m_faultOutputManager[i]->setInputParam(seissolInstance.meshReader());
      m_faultOutputManager[i]->setLtsData(
          &m_ltsTree, &m_lts, &m_ltsLut, m_dynRupTree[i], m_dynRup[i].get());
      m_faultOutputManager[i]->setBackupTimeStamp(backupTimeStamp);
      m_faultOutputManager[i]->init();
    }
  }
}


void seissol::initializer::MemoryManager::initFrictionData() {
  for(int i=0; i < MULTIPLE_SIMULATIONS; i++){

  if (m_seissolParams->drParameters[i]->isDynamicRuptureEnabled) {

    m_DRInitializer[i]->initializeFault(m_dynRup[i].get(), m_dynRupTree[i]);

#ifdef ACL_DEVICE
    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(m_FrictionLawDevice.get())) {

      const auto mask = seissol::initializer::LayerMask(Ghost);

      impl->allocateAuxiliaryMemory();
    }
#endif // ACL_DEVICE
  }
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

