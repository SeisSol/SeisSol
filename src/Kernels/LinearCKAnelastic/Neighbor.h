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

#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"

namespace seissol::kernels::solver::linearckanelastic {
class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData<Config>& global) override;

  void computeNeighborsIntegral(
      LTS::Ref<Config>& data,
      const std::array<real*, Cell::NumFaces>& timeIntegrated,
      const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) override;

  void computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] std::pair<PerformanceEstimate, PerformanceEstimate> metrics(
      const std::array<FaceType, Cell::NumFaces>& faceTypes,
      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
      const std::array<CellDRMapping<Config>, Cell::NumFaces>& cellDrMapping) const override;

  protected:
  kernel::neighborFluxExt<Config> nfKrnlPrototype_;
  kernel::neighbor<Config> nKrnlPrototype_;
  dynamicRupture::kernel::nodalFlux<Config> drKrnlPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_neighborFluxExt<Config> deviceNfKrnlPrototype_;
  kernel::gpu_neighbor<Config> deviceNKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux<Config> deviceDrKrnlPrototype_;
#endif
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_NEIGHBOR_H_
