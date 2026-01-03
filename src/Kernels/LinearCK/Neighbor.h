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
#include "GeneratedCode/kernel.h"
#include "Kernels/Neighbor.h"
#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels::solver::linearck {

class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(LTS::Ref& data,
                                const CellDRMapping (&cellDrMapping)[4],
                                real* timeIntegrated[4],
                                real* faceNeighborsPrefetch[4]) override;

  void computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
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
  kernel::neighboringFlux nfKrnlPrototype_;
  dynamicRupture::kernel::nodalFlux drKrnlPrototype_;

#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_NEIGHBOR_H_
