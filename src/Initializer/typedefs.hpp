/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2020, SeisSol Group
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

#include "BasicTypedefs.hpp"
#include <Initializer/preProcessorMacros.fpp>
#include <Kernels/equations.hpp>
#include "Equations/datastructures.hpp"
#include <generated_code/tensor.h>

#include <cstddef>

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
  FaceType faceTypes[4];

  // mapping of the neighboring elements to the references element in relation to this element
  int faceRelations[4][2];

  // ids of the face neighbors
  unsigned int faceNeighborIds[4];

  // LTS setup
  unsigned short ltsSetup;

  // unique global id of the time cluster
  unsigned int clusterId;
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
   * Addresses of the global change of basis matrices (multiplied by the inverse diagonal mass matrix):
   * 
   *    0: \f$ M^{-1} R^1 \f$
   *    1: \f$ M^{-1} R^2 \f$
   *    2: \f$ M^{-1} R^3 \f$
   *    3: \f$ M^{-1} R^4 \f$
   **/
  seissol::tensor::rDivM::Container<real const*> changeOfBasisMatrices;
  
  /**
   * Addresses of the transposed global change of basis matrices left-multiplied with the local flux matrix:
   * 
   *    0: \f$ F^- ( R^1 )^T \f$
   *    1: \f$ F^- ( R^2 )^T \f$
   *    2: \f$ F^- ( R^3 )^T \f$
   *    3: \f$ F^- ( R^4 )^T \f$
   **/
  seissol::tensor::fMrT::Container<real const*> localChangeOfBasisMatricesTransposed;
  
  /**
   * Addresses of the transposed global change of basis matrices:
   * 
   *    0: \f$ ( R^1 )^T \f$
   *    1: \f$ ( R^2 )^T \f$
   *    2: \f$ ( R^3 )^T \f$
   *    3: \f$ ( R^4 )^T \f$
   **/
  seissol::tensor::rT::Container<real const*> neighbourChangeOfBasisMatricesTransposed;
  
  /**
   * Addresses of the global flux matrices:
   * 
   *    0: \f$ F^{+,1} \f$
   *    1: \f$ F^{+,2} \f$
   *    2: \f$ F^{+,3} \f$
   **/
  seissol::tensor::fP::Container<real const*> neighbourFluxMatrices;

  /** 
   * Addresses of the global stiffness matrices (multiplied by the inverse diagonal mass matrix):
   *
   *    0:  \f$ M^{-1} K^\xi \f$
   *    1:  \f$ M^{-1} K^\eta \f$
   *    2:  \f$ M^{-1} K^\zeta f$
   *
   *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks (except for the additional flux matrix).
   **/ 
  seissol::tensor::kDivM::Container<real const*> stiffnessMatrices;

  /** 
   * Addresses of the transposed global stiffness matrices (multiplied by the inverse diagonal mass matrix):
   *
   *    0:  \f$ M^{-1} ( K^\xi )^T \f$
   *    1:  \f$ M^{-1} ( K^\eta )^T \f$
   *    2:  \f$ M^{-1} ( K^\zeta )^T \f$
   *
   *   Remark: The ordering of the pointers is identical to the ordering of the memory chunks (except for the additional flux matrix).
   **/ 
  seissol::tensor::kDivMT::Container<real const*> stiffnessMatricesTransposed;

  /**
   * Address of the (thread-local) local time stepping integration buffers used in the neighbor integral computation
   **/
  real *integrationBufferLTS{nullptr};
  
   /** 
   * Addresses of the global nodal flux matrices
   *
   *    0:  \f$ P^{+,1} \f$
   *    1:  \f$ P^{-,1,1} \f$
   *    2:  \f$ P^{-,1,2} \f$
   *    3 : \f$ P^{-,1,3} \f$
   *    4:  \f$ P^{+,2} \f$
   *    5:  \f$ P^{-,2,1} \f$
   *    6:  \f$ P^{-,2,2} \f$
   *    7 : \f$ P^{-,2,3} \f$
   *    [..]
   *    15: \f$ P^{-,4,3} \f$
   **/ 
  seissol::tensor::V3mTo2nTWDivM::Container<real const*> nodalFluxMatrices;

  seissol::nodal::tensor::V3mTo2nFace::Container<real const*> V3mTo2nFace;
  seissol::tensor::project2nFaceTo3m::Container<real const*> project2nFaceTo3m;

  /** 
   * Addresses of the global face to nodal matrices
   *
   *    0:  \f$ N^{+,1} \f$
   *    1:  \f$ N^{-,1,1} \f$
   *    2:  \f$ N^{-,1,2} \f$
   *    3 : \f$ N^{-,1,3} \f$
   *    4:  \f$ N^{+,2} \f$
   *    5:  \f$ N^{-,2,1} \f$
   *    6:  \f$ N^{-,2,2} \f$
   *    7 : \f$ N^{-,2,3} \f$
   *    [..]
   *    15: \f$ N^{-,4,3} \f$
   **/ 

 
  seissol::tensor::V3mTo2n::Container<real const*> faceToNodalMatrices;

  //! Modal basis to quadrature points
  real* evalAtQPMatrix{nullptr};

  //! Project function evaluated at quadrature points to modal basis
  real* projectQPMatrix{nullptr};
  
  //! Switch to nodal for plasticity
  real* vandermondeMatrix{nullptr};
  real* vandermondeMatrixInverse{nullptr};

  // A vector of ones. Note: It is only relevant for GPU computing.
  // It allows us to allocate this vector only once in the GPU memory
  real* replicateStresses{nullptr};
};

struct CompoundGlobalData {
  GlobalData* onHost{nullptr};
  GlobalData* onDevice{nullptr};
};

// data for the cell local integration
struct LocalIntegrationData {
  // star matrices
  real starMatrices[3][seissol::tensor::star::size(0)];

  // flux solver for element local contribution
  real nApNm1[4][seissol::tensor::AplusT::size()];

  // equation-specific data
  //TODO(Lukas/Sebastian):
  //Get rid of ifdefs
#if defined USE_ANISOTROPIC
  seissol::model::AnisotropicLocalData specific;
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  seissol::model::ViscoElasticLocalData specific;
#elif defined USE_ELASTIC
  seissol::model::ElasticLocalData specific;
#endif
};

// data for the neighboring boundary integration
struct NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  real nAmNm1[4][seissol::tensor::AminusT::size()];

  // equation-specific data
  //TODO(Lukas/Sebastian):
  //Get rid of ifdefs
#if defined USE_ANISOTROPIC
  seissol::model::AnisotropicNeighborData specific;
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  seissol::model::ViscoElasticNeighborData specific;
#elif defined USE_ELASTIC
  seissol::model::ElasticNeighborData specific;
#endif
};

// material constants per cell
struct CellMaterialData {
  //TODO(Lukas/Sebastian):
  //Get rid of ifdefs
#if defined USE_ANISOTROPIC
  seissol::model::AnisotropicMaterial local;
  seissol::model::AnisotropicMaterial neighbor[4];
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
  seissol::model::ViscoElasticMaterial local;
  seissol::model::ViscoElasticMaterial neighbor[4];
#elif defined USE_ELASTIC
  seissol::model::ElasticMaterial local;
  seissol::model::ElasticMaterial neighbor[4];
#else
  static_assert(false, "No Compiler flag for the material behavior has been given. Current implementation allows: USE_ANISOTROPIC, USE_ISOTROPIC, USE_VISCOELASTIC, USE_VISCOELASTIC2");
#endif
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

struct DRFaceInformation {
  unsigned meshFace;
  unsigned plusSide;
  unsigned minusSide;
  unsigned faceRelation;
};

struct DRGodunovData {
  real TinvT[seissol::tensor::TinvT::size()];
};

struct CellDRMapping {
  unsigned side;
  unsigned faceRelation;
  real* godunov;
  real* fluxSolver;
};

struct CellBoundaryMapping {
  real* nodes;
  real* TData;
  real* TinvData;
  real* easiBoundaryConstant;
  real* easiBoundaryMap;
};

struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  real nodes[seissol::nodal::tensor::nodes2D::Shape[0] * 3];
  real TData[seissol::tensor::T::size()];
  real TinvData[seissol::tensor::Tinv::size()];
  real easiBoundaryConstant[seissol::tensor::easiBoundaryConstant::size()];
  real easiBoundaryMap[seissol::tensor::easiBoundaryMap::size()];
};

/*
 * \class MemoryProperties
 *
 * \brief An auxiliary data structure for a policy-based design
 *
 * Attributes are initialized with CPU memory properties by default.
 * See, an example of a policy-based design in GlobalData.cpp
 * */
struct MemoryProperties {
  size_t alignment{ALIGNMENT};
  size_t pagesizeHeap{PAGESIZE_HEAP};
  size_t pagesizeStack{PAGESIZE_STACK};
};

namespace seissol {
struct GravitationSetup {
  double acceleration = 9.81; // m/s
};
} // namespace seissol

struct TravellingWaveParameters {
  std::array<double, 3> origin;
  std::array<double, 3> kVec;
  std::vector<int> varField;
  std::vector<std::complex<double>> ampField;
};

#endif
