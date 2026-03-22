// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_NEIGHBOR_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/Kernel.h"
#include "Parallel/Runtime/Stream.h"

namespace seissol::kernels {

class NeighborKernel : public Kernel {
  public:
  ~NeighborKernel() override = default;

  /**
   * @brief Compute the flux contribution of neighboring cells onto the cell, given time-integrated
   * neighbor cell DoFs.
   *
   * This step equals the second part of the ADER-DG "corrector"; meaning that we add the flux
   * contributions from the neighboring cells as well as the dynamic rupture faces. We assume that
   * we have the neighboring values already present in a time-integrated manner. We either take them
   * from the neighboring cells directly (they are usually precomputed), or need to have them
   * pre-computed via the TimeKernel (mostly in the case of LTS, if our neighbor is a cell from a
   * larger time cluster).
   *
   * @param data Cell data reference object (contains references to all stored data arrays for that
   * cell)
   * @param timeIntegrated The time-integrated DoFs of neighboring cells (usually ONLY for regular
   * faces; DR is handled via the data object normally)
   * @param faceNeighborsPrefetch The current time step width
   */
  virtual void
      computeNeighborsIntegral(LTS::Ref& data,
                               const std::array<real*, Cell::NumFaces>& timeIntegrated,
                               const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) = 0;

  virtual void
      computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsNeighborsIntegral(
      const std::array<FaceType, Cell::NumFaces>& faceTypes,
      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
      const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping,
      std::uint64_t& nonZeroFlops,
      std::uint64_t& hardwareFlops,
      std::uint64_t& drNonZeroFlops,
      std::uint64_t& drHardwareFlops) = 0;

  virtual std::uint64_t bytesNeighborsIntegral() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_NEIGHBOR_H_
