// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
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
#include "Kernels/NeighborBase.h"
#include "Parallel/Runtime/Stream.h"

namespace seissol::kernels {

class Neighbor : public NeighborBase {
  public:
  void setHostGlobalData(const GlobalData* global);
  void setGlobalData(const CompoundGlobalData& global);

  void computeNeighborsIntegral(NeighborData& data,
                                const CellDRMapping (&cellDrMapping)[4],
                                real* timeIntegrated[4],
                                real* faceNeighborsPrefetch[4]);

  void computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsNeighborsIntegral(const FaceType faceTypes[4],
                              const int neighboringIndices[4][2],
                              const CellDRMapping (&cellDrMapping)[4],
                              unsigned int& nonZeroFlops,
                              unsigned int& hardwareFlops,
                              long long& drNonZeroFlops,
                              long long& drHardwareFlops);

  unsigned bytesNeighborsIntegral();
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_NEIGHBOR_H_
