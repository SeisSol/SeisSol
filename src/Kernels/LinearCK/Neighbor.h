// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_

#include "Common/Constants.h"
#include "Common/Real.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"
#include "Monitoring/Metric.h"
#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::linearck {

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
  kernel::neighboringFlux<Cfg> nfKrnlPrototype_;
  dynamicRupture::kernel::nodalFlux<Cfg> drKrnlPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux<Cfg> deviceNfKrnlPrototype_;
  dynamicRupture::kernel::gpu_nodalFlux<Cfg> deviceDrKrnlPrototype_;
  device::DeviceInstance& device_ = device::DeviceInstance::instance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_
