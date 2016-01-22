/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 * Typedefs for the implementation.
 **/

#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <Initializer/preProcessorMacros.fpp>
#include <Kernels/precision.hpp>
#include <Kernels/equations.hpp>
#include <Model/datastructures.hpp>
#include <generated_code/sizes.h>

#include <cstddef>

enum Layer {
  ghost,
  copy,
  interior
};

enum mpiTag {
  localIntegrationData = 0,
  neighboringIntegrationData = 1,
  timeData = 2
};

enum TimeClustering {
  // global time stepping
  single    = 0,
  // offline clustering computed in pre-processing
  offline   = 1,
  // online clustering resulting in a multi-rate scheme
  multiRate = 2,
  // online clustering aiming at LTS for slithers only
  slithers  = 3
};

// face types
enum faceType {
  // regular: inside the computational domain
  regular = 0,

  // free surface boundary
  freeSurface = 1,
  
  // dynamic rupture boundary
  dynamicRupture = 3,

  // absorbing/outflow boundary
  outflow = 5,

  // periodic boundary
  periodic = 6
};

// cross-cluster time stepping information
struct TimeStepping {
  /*
   * Number of lts clusters prensent throughout the entire domain.
   */
  unsigned int numberOfGlobalClusters;

  /*
   * Time step rate the cluster is smaller than the next one.
   */
  unsigned int *globalTimeStepRates;

  /*
   * Time step widths of the global time clusters according to CFL.
   */
  double *globalCflTimeStepWidths;

  /*
   * Time all clusters are aiming to reach.
   */
  double synchronizationTime;

  /*
   * Number of local clusters present on this rank.
   */
  unsigned int numberOfLocalClusters;

  /*
   * Ids of the local clusters with respect to global ordering.
   */
  unsigned int *clusterIds;
};

// cell local information
struct CellLocalInformation {
  // types of the faces
  enum faceType faceTypes[4];

  // mapping of the neighboring elements to the references element in relation to this element
  int faceRelations[4][2];

  // ids of the face neighbors
  unsigned int faceNeighborIds[4];

  // LTS setup
  unsigned short ltsSetup;

  // unique global id of the time cluster
  unsigned int clusterId;

  // maximum cfl time step width of this cell
  double timeStepWidth;
};

struct MeshStructure {
  /*
   * Number of regions in the ghost and copy layer.
   * This is equivalent to the number of ranks in a MPI setting.
   */
  unsigned int numberOfRegions;

  /*
   * Region-specific neighboring clusters
   * [0]: rank
   * [1]: global time cluster id
   */
  int (*neighboringClusters)[2];

  /*
   * Total number of ghost cells.
   */
  unsigned int numberOfGhostCells;

  /*
   * Number of ghost cells in each region of the ghost layer.
   */
  unsigned int *numberOfGhostRegionCells;

  /*
   * Number of cells with derivatives in each region of the ghost layer.
   */
  unsigned int *numberOfGhostRegionDerivatives;

  /*
   * Pointers to the memory chunks of the ghost regions.
   */
  real** ghostRegions;

  /*
   * Sizes of the ghost regions (in reals).
   */
  unsigned int *ghostRegionSizes;

  /*
   * Total number of copy cells.
   */
  unsigned int numberOfCopyCells;

  /*
   * Number of copy cells in each region of the copy layer.
   */
  unsigned int *numberOfCopyRegionCells;

  /*
   * Number of cells with communicating derivatives in each region of the ghost layer.
   */
  unsigned int *numberOfCommunicatedCopyRegionDerivatives;

  /*
   * Pointers to the memory chunks of the copy regions.
   *   Remark: For the cells in the copy layer more information will be stored (in general).
   *           The pointers only point to communcation related chunks.
   */
  real** copyRegions;

  /*
   * Sizes of the copy regions (in reals).
   */
  unsigned int *copyRegionSizes;


  /*
   * Total number of interior cells without MPI-face-neighbors.
   */
  unsigned int numberOfInteriorCells;

  /*
   * Message identifiers for the sends.
   */
  int *sendIdentifiers;

  /*
   * Message identifiers for the receives.
   */
  int *receiveIdentifiers;

#ifdef USE_MPI
  /*
   * MPI send requests.
   */
  MPI_Request *sendRequests;

  /*
   * MPI receive requests.
   */
  MPI_Request *receiveRequests;
#endif

};

struct GlobalData {
 /** 
   * Addresses of the global flux matrices (multiplied by the inverse diagonal mass matrix):
   *
   *    0:  \f$ M^{-1} F^{-, 1} \f$
   *    1 : \f$ M^{-1} F^{-, 2} \f$
   *    2:  \f$ M^{-1} F^{-, 3} \f$
   *    3 : \f$ M^{-1} F^{-, 4} \f$
   *    4:  \f$ M^{-1} F^+{+, 1, 1, 1} \f$
   *    5:  \f$ M^{-1} F^+{+, 1, 1, 2} \f$
   *    6:  \f$ M^{-1} F^+{+, 1, 1, 3} \f$
   *    7:  \f$ M^{-1} F^+{+, 1, 2, 1} \f$
   *    8:  \f$ M^{-1} F^+{+, 1, 1, 2} \f$
   *    9:  \f$ M^{-1} F^+{+, 1, 1, 3} \f$
   *    [..]
   *    51: \f$ M^{-1} F^+{+, 4, 4, 3} \f$
   *    52: \f$ N_{k,i} A_k^+ N_{k,i}^{-1}\f$ or \f$ N_{k,i} A_{k(i)}^- N_{k,i}^{-1} \f$
   *
   *   Remark: The ordering of the pointers is given as above, however the chunks in memory are allowed to have a different orderning as given in the XML-file.
   **/ 
  real *fluxMatrices[52];

  /** 
   * Addresses of the global stiffness matrices (multiplied by the inverse diagonal mass matrix):
   *
   *    0:  \f$ M^{-1} K^\xi \f$
   *    1:  \f$ M^{-1} K^\eta \f$
   *    2:  \f$ M^{-1} K^\zeta f$
   *
   *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks (except for the additional flux matrix).
   **/ 
  real *stiffnessMatrices[3];

  /** 
   * Addresses of the transposed global stiffness matrices (multiplied by the inverse diagonal mass matrix):
   *
   *    0:  \f$ M^{-1} ( K^\xi )^T \f$
   *    1:  \f$ M^{-1} ( K^\eta )^T \f$
   *    2:  \f$ M^{-1} ( K^\zeta )^T \f$
   *
   *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks (except for the additional flux matrix).
   **/ 
  real *stiffnessMatricesTransposed[3];
  
  /**
   * Address of the global inverse mass matrix
   **/
  real *inverseMassMatrix;

  /**
   * Address of the (thread-local) local time stepping integration buffers used in the neighbor integral computation
   **/
  real *integrationBufferLTS;
};

// data for the cell local integration
struct LocalIntegrationData {
  // star matrices
  real starMatrices[3][seissol::model::AstarT::reals];

  // flux solver for element local contribution
  real nApNm1[4][seissol::model::AplusT::reals];

  // equation-specific data
  seissol::model::LocalData specific;
};

// data for the neighboring boundary integration
struct NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  real nAmNm1[4][seissol::model::AminusT::reals];

  // equation-specific data
  seissol::model::NeighborData specific;
};

// material constants per cell
struct CellMaterialData {
  seissol::model::Material local;
  seissol::model::Material neighbor[4];
};

// plasticity information per cell
struct PlasticityData {
  // initial loading (stress tensor)
  real initialLoading[6][NUMBER_OF_BASIS_FUNCTIONS];
};

/**
 * Cell local data (material dependent).
 **/
struct CellData {
  // local integration data
  struct LocalIntegrationData       *localIntegration;
  // neighboring integration data
  struct NeighboringIntegrationData *neighboringIntegration;
  // local and neighbor material data
  CellMaterialData                  *material;
  // Plasticity
  PlasticityData                    *plasticity;
};

/**
 * Internal state of the wave propagation component.
 *
 * Splits into:
 *
 *  * ghostLayer:    time buffers and/or derivatives of neighboring ranks required for cell updates in the computational domain.
 *                   stride-1 memory access per rank is enforced.
 *
 *  * copyLayer:     time buffers and/or derivatives required for updates on neighboring ranks and in this ranks computational domain.
 *                   stride-1 memory access for communication is enforced, this includes duplicated cells if the information is required for more than one neighbor.
 *
 *  * interiorTime:  time buffers and/or derivatives required for updates in the computational domain only.
 *
 *  * buffers /      pointers for every cells, ordering: 1) layer (ghost, copy, interior) 2) time cluster.
 *    derivatives
 *    face neighbors
 *
 *  * dofs:          modal degrees of freedom for all cells in the computational domain (copy + interior).
 *
 **/
struct InternalState {
#ifdef USE_MPI
  /*
   * Ghost layer: Buffers and derivatives
   */
  real (*ghostLayer);

  /*
   * Copy layer: Buffers and derivatives
   */
  real (*copyLayer);
#endif

  /*
   * Time Interior: Buffers and derivatives
   */
  real (*interiorTime);

  /*
   * Pointers to time buffers (or NULL if not present).
   *   One pointer per cell (including duplicate cells in the copy layer).
   *   Covers ghost layer, copy layer and interior.
   *   Sorting analogue to mesh (cross-cluster cell ids are valid):
   *     1) Time cluster id
   *     2) Layer: Ghost, copy, interior
   */
  real **buffers;

  /*
   * Pointers to time derivatives (or NULL if not present).
   *   One pointer per cell (including duplicate cells in the copy layer).
   *   Covers ghost layer, copy layer and interior.
   *   Sorting analogue to mesh (cross-cluster cell ids are valid):
   *     1) Time cluster id
   *     2) Layer: Ghost, copy, interior
   */
  real **derivatives;

  /*
   * Pointers to the either the time buffers or time derivatives of the face neighbors.
   *   One pointer per cell (including duplicates in the copy layer).
   *   Covers copy layer and interior (excludes pointers for the ghost layer).
   *   Sorting analogue to mesh:
   *     1) Time cluster id
   *     2) Layer: Copy, interior
   */
  real *(*faceNeighbors)[4];

  /*
   * Regular degrees of freedoom in copy and interior: modal coefficients.
   *   Covers copy layer and interior (excludes ghost layer).
   */
  real (*dofs)[NUMBER_OF_ALIGNED_DOFS];

  // plastic strain
  real (*pstrain)[7];
};

/**
 * Structure of the cells.
 **/
struct Cells {
#ifdef USE_MPI
  /*
   * Regular degrees of freedom in the copy layer: modal coefficients.
   */
  real (*copyDofs)[NUMBER_OF_ALIGNED_DOFS];

  /*
   * Pointers to time buffers in the copy layer (or NULL if not used).
   */
  real **copyBuffers;

  /*
   * Pointers to derivatives int the copy layer (or NULL if not used).
   */
  real **copyDerivatives;

  /*
   * Pointers to the either the time buffers or time derivatives of the face neighbors in the copy layer.
   */
  real *(*copyFaceNeighbors)[4];

  /** Pointer to copy layer plastic strain */
  real (*copyPstrain)[7];
#endif

  /*
   * Regular degrees of freedom in the interior: modal coefficients.
   */
  real (*interiorDofs)[NUMBER_OF_ALIGNED_DOFS];

  /*
   * Pointers to time buffers in the interior (or NULL if not used).
   */
  real **interiorBuffers;

  /*
   * Pointers to derivatives int the interior (or NULL if not used).
   */
  real **interiorDerivatives;

  /*
   * Pointers to the either the time buffers or time derivatives of the face neighbors in the interior.
   */
  real *(*interiorFaceNeighbors)[4];

  /** Pointer to interior plastic strain */
  real (*interiorPstrain)[7];
};

/** A piecewise linear function.
 * 
 *  Say t \in I_j, then
 *    f(t) = m_j * t + n_j,
 *  where I_j is the half-open interval [t_o + j*dt, t_o + (j+1)*dt).
 *  j runs through 0,...,n-1.
 **/
struct PiecewiseLinearFunction1D {
   /** slopes[i] = m_i */
  real* slopes;

  /** intercepts[i] = n_i */
  real* intercepts;
  
  /** numberOfPieces = n */
  unsigned numberOfPieces;
  
  /** onsetTime = t_o */
  real onsetTime;
  
  /** samplingInterval = dt */
  real samplingInterval;
  
  PiecewiseLinearFunction1D() : slopes(NULL), intercepts(NULL), numberOfPieces(0) {}
  ~PiecewiseLinearFunction1D() { delete[] slopes; delete[] intercepts; numberOfPieces = 0; }
};

#endif
