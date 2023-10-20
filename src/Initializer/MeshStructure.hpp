// Copyright (c) 2015-2020 SeisSol Group
// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MESHSTRUCTURE_20231020_HPP
#define MESHSTRUCTURE_20231020_HPP

#include <Kernels/precision.hpp>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace seissol {

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
  unsigned int* numberOfGhostRegionCells;

  /*
   * Number of cells with derivatives in each region of the ghost layer.
   */
  unsigned int* numberOfGhostRegionDerivatives;

  /*
   * Pointers to the memory chunks of the ghost regions.
   */
  real** ghostRegions;

  /*
   * Sizes of the ghost regions (in reals).
   */
  unsigned int* ghostRegionSizes;

  /*
   * Total number of copy cells.
   */
  unsigned int numberOfCopyCells;

  /*
   * Number of copy cells in each region of the copy layer.
   */
  unsigned int* numberOfCopyRegionCells;

  /*
   * Number of cells with communicating derivatives in each region of the ghost layer.
   */
  unsigned int* numberOfCommunicatedCopyRegionDerivatives;

  /*
   * Pointers to the memory chunks of the copy regions.
   *   Remark: For the cells in the copy layer more information will be stored (in general).
   *           The pointers only point to communcation related chunks.
   */
  real** copyRegions;

  /*
   * Sizes of the copy regions (in reals).
   */
  unsigned int* copyRegionSizes;

  /*
   * Total number of interior cells without MPI-face-neighbors.
   */
  unsigned int numberOfInteriorCells;

  /*
   * Message identifiers for the sends.
   */
  int* sendIdentifiers;

  /*
   * Message identifiers for the receives.
   */
  int* receiveIdentifiers;

#ifdef USE_MPI
  /*
   * MPI send requests.
   */
  MPI_Request* sendRequests;

  /*
   * MPI receive requests.
   */
  MPI_Request* receiveRequests;
#endif
};

} // namespace seissol

#endif // MESHSTRUCTURE_20231020_HPP
