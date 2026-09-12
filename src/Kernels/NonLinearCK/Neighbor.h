// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_NEIGHBOR_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Neighbor.h"
#include "Monitoring/Metric.h"

#include <array>
#include <cstdint>
#include <utility>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels::solver::nonlinearck {

/// The far side of each face.
///
/// The neighbour hands over the same two integrals every cell produces, the
/// state and the stress, so this kernel evaluates its half of the Rusanov flux
/// without ever reading the neighbour's material. Only one number about the
/// other cell is needed, its largest wave speed, and it arrives through the
/// neighbour data rather than through the material.
///
/// Dynamic rupture faces are not handled. The material declares that it does
/// not support them, because the traction a friction law needs is derived from
/// the strain here rather than stored.
class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(
      LTS::Ref& data,
      const std::array<real*, Cell::NumFaces>& timeIntegrated,
      const std::array<real*, Cell::NumFaces>& faceNeighborsPrefetch) override;

  void computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] std::pair<PerformanceEstimate, PerformanceEstimate>
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
              const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
              const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const override;

  protected:
  kernel::projectToFace projectToFace_;
  kernel::projectNeighborToFace projectNeighborToFace_;
  kernel::damageRusanov rusanov_;
  kernel::faceIntegral faceIntegral_;

#ifdef ACL_DEVICE
  kernel::gpu_projectToFace deviceProjectToFace_;
  kernel::gpu_projectNeighborToFace deviceProjectNeighborToFace_;
  kernel::gpu_damageRusanov deviceRusanov_;
  kernel::gpu_faceIntegral deviceFaceIntegral_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::nonlinearck

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_NEIGHBOR_H_
