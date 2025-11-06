// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_

#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"

namespace seissol::kernels::solver::linearckanelastic {
class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(LTS::Ref& data,
                                const CellDRMapping (&cellDrMapping)[4],
                                real* timeIntegrated[4],
                                real* faceNeighborsPrefetch[4]) override;

  void computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsNeighborsIntegral(
      const std::array<FaceType, Cell::NumFaces>& faceTypes,
      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
      const CellDRMapping (&cellDrMapping)[4],
      std::uint64_t& nonZeroFlops,
      std::uint64_t& hardwareFlops,
      std::uint64_t& drNonZeroFlops,
      std::uint64_t& drHardwareFlops) override;

  std::uint64_t bytesNeighborsIntegral() override;

  protected:
  kernel::neighborFluxExt m_nfKrnlPrototype;
  kernel::neighbor m_nKrnlPrototype;
  dynamicRupture::kernel::nodalFlux m_drKrnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt deviceNfKrnlPrototype;
  kernel::gpu_neighbor deviceNKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_
