// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBORBASE_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBORBASE_H_

#include "GeneratedCode/kernel.h"
#include <Kernels/Neighbor.h>
#include <Memory/GlobalData.h>

namespace seissol::kernels::solver::linearckanelastic {

template <typename Cfg>
class Neighbor : public NeighborKernel<Cfg> {
  public:
  using real = Real<Cfg>;

  void setGlobalData(const GlobalData& global) override;

  void computeNeighborsIntegral(LTS::Ref<Cfg>& data,
                                const CellDRMapping<Cfg> (&cellDrMapping)[4],
                                real* timeIntegrated[4],
                                real* faceNeighbors_prefetch[4]) override;

  void computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsNeighborsIntegral(
      const std::array<FaceType, Cell::NumFaces>& faceTypes,
      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
      const CellDRMapping<Cfg> (&cellDrMapping)[4],
      std::uint64_t& nonZeroFlops,
      std::uint64_t& hardwareFlops,
      std::uint64_t& drNonZeroFlops,
      std::uint64_t& drHardwareFlops) override;

  std::uint64_t bytesNeighborsIntegral() override;

  protected:
  kernel::neighborFluxExt<Cfg> m_nfKrnlPrototype;
  kernel::neighbor<Cfg> m_nKrnlPrototype;
  dynamicRupture::kernel::nodalFlux<Cfg> m_drKrnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt<Cfg> deviceNfKrnlPrototype;
  kernel::gpu_neighbor<Cfg> deviceNKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> deviceDrKrnlPrototype;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBORBASE_H_
