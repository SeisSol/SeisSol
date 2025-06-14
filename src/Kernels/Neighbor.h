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
#include "Parallel/Runtime/Stream.h"
#include <Kernels/Kernel.h>

namespace seissol::kernels {

class NeighborKernel : public Kernel {
  public:
  ~NeighborKernel() override = default;

  virtual void computeNeighborsIntegral(NeighborData& data,
                                        const CellDRMapping (&cellDrMapping)[4],
                                        real* timeIntegrated[4],
                                        real* faceNeighborsPrefetch[4]) = 0;

  virtual void
      computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsNeighborsIntegral(const FaceType faceTypes[4],
                                      const int neighboringIndices[4][2],
                                      const CellDRMapping (&cellDrMapping)[4],
                                      unsigned int& nonZeroFlops,
                                      unsigned int& hardwareFlops,
                                      long long& drNonZeroFlops,
                                      long long& drHardwareFlops) = 0;

  virtual unsigned bytesNeighborsIntegral() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_NEIGHBOR_H_
