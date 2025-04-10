// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_

#include "generated_code/kernel.h"
#include <Kernels/Neighbor.h>

namespace seissol::kernels::solver::linearckanelastic {
class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(NeighborData& data,
                                const CellDRMapping (&cellDrMapping)[4],
                                real* timeIntegrated[4],
                                real* faceNeighborsPrefetch[4]) override;

  void computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsNeighborsIntegral(const FaceType faceTypes[4],
                              const int neighboringIndices[4][2],
                              const CellDRMapping (&cellDrMapping)[4],
                              unsigned int& nonZeroFlops,
                              unsigned int& hardwareFlops,
                              long long& drNonZeroFlops,
                              long long& drHardwareFlops) override;

  unsigned bytesNeighborsIntegral() override;

  protected:
  kernel::neighbourFluxExt m_nfKrnlPrototype;
  kernel::neighbour m_nKrnlPrototype;
  dynamicRupture::kernel::nodalFlux m_drKrnlPrototype;
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_
