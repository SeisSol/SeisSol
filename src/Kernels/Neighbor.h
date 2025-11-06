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

  virtual void computeNeighborsIntegral(LTS::Ref& data,
                                        const CellDRMapping (&cellDrMapping)[4],
                                        real* timeIntegrated[4],
                                        real* faceNeighborsPrefetch[4]) = 0;

  virtual void
      computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsNeighborsIntegral(
      const std::array<FaceType, Cell::NumFaces>& faceTypes,
      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
      const CellDRMapping (&cellDrMapping)[4],
      std::uint64_t& nonZeroFlops,
      std::uint64_t& hardwareFlops,
      std::uint64_t& drNonZeroFlops,
      std::uint64_t& drHardwareFlops) = 0;

  virtual std::uint64_t bytesNeighborsIntegral() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_NEIGHBOR_H_
