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

#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"

namespace seissol::kernels::solver::linearckanelastic {
template <typename Cfg>
class Neighbor : public NeighborKernel<Cfg> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  void setGlobalData(const CompoundGlobalData<Cfg>& global) override;

  void computeNeighborsIntegral(
      LTS::Ref<Cfg>& data,
      const std::array<real*, Cell::NumFaces>& timeIntegrated,
      const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) override;

  void computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] std::pair<PerformanceEstimate, PerformanceEstimate>
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
              const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
              const std::array<CellDRMapping<Cfg>, Cell::NumFaces>& cellDrMapping) const override;

  protected:
  kernel::neighborFluxExt<Cfg> nfKrnlPrototype_;
  kernel::neighbor<Cfg> nKrnlPrototype_;
  dynamicRupture::kernel::nodalFlux<Cfg> drKrnlPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt<Cfg> deviceNfKrnlPrototype_;
  kernel::gpu_neighbor<Cfg> deviceNKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> deviceDrKrnlPrototype_;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_
