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
/// A rupture face is the one kind that contributes nothing through this pair:
/// its state does not cross it, a friction law decides what does. What is left
/// for this kernel is to lift that decision back onto the cell, which is the
/// same nodal flux the linear solver applies -- what differs is the flux solver
/// it is applied with, and that was built at setup.
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
  kernel::damageLocalFlux localFlux_;
  kernel::damageNeighborFlux neighborFlux_;
  dynamicRupture::kernel::nodalFlux drFlux_;

#ifdef ACL_DEVICE
  kernel::gpu_damageLocalFlux deviceLocalFlux_;
  kernel::gpu_damageNeighborFlux deviceNeighborFlux_;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrFlux_;
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::nonlinearck

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_NEIGHBOR_H_
