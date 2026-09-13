// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Neighbor.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Common/Offset.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/Typedefs.h"
#include "Kernels/NonLinearCK/Solver.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <utils/logger.h>
#include <yateto.h>

namespace seissol::kernels::solver::nonlinearck {

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
  localFlux_.bindGlobals(*global.onHost);
  neighborFlux_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceLocalFlux_.bindGlobals(*global.onDevice);
  deviceNeighborFlux_.bindGlobals(*global.onDevice);
#endif
}

void Neighbor::computeNeighborsIntegral(
    LTS::Ref& data,
    const std::array<real*, Cell::NumFaces>& timeIntegrated,
    const std::array<real*, Cell::NumFaces>& /*faceNeighborsPrefetch*/) {
  const auto& info = data.get<LTS::CellInformation>();

  // The cell's own integrals: the ones its predictor wrote and its neighbours
  // read from it. That it has them is not luck -- the solver declares it needs
  // them and the LTS setup gives every cell the buffer for it, including one
  // whose faces are all boundaries.
  const real* own = data.get<LTS::StepIntegrals>();
  assert(own != nullptr && Solver::RequiresOwnIntegrals);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const auto faceType = info.faceTypes[face];
    const bool hasNeighbor = faceType == FaceType::Regular;

    // A face without a neighbour carries its ghost rule in its own pair of
    // matrices, folded in at setup. What is left of it here is that there is
    // no far state: its own stands in, which is also the speed the face is
    // scaled with.
    const real* other = hasNeighbor ? timeIntegrated[face] : own;

    // Each half assembles the pair it needs and reads the face's speed out of
    // the two tensors itself. Nothing about a face crosses a kernel boundary,
    // which is what lets the batched path run on the same tables the linear
    // solver records.
    kernel::damageLocalFlux local = localFlux_;
    local.Q = data.get<LTS::Dofs>();
    local.I = own;
    local.INeighbor = other;
    local.fluxConstant = data.get<LTS::LocalIntegration>().nApNm1[face];
    local.fluxDissipation = data.get<LTS::NeighboringIntegration>().nAmNm1[face];
    local.execute(face);

    if (!hasNeighbor) {
      continue;
    }

    kernel::damageNeighborFlux neighbor = neighborFlux_;
    neighbor.Q = data.get<LTS::Dofs>();
    neighbor.I = own;
    neighbor.INeighbor = other;
    neighbor.fluxConstant = data.get<LTS::LocalIntegration>().nApNm1[face];
    neighbor.fluxDissipation = data.get<LTS::NeighboringIntegration>().nAmNm1[face];
    neighbor.execute(info.faceRelations[face][1], info.faceRelations[face][0], face);
  }
}

void Neighbor::computeBatchedNeighborsIntegral(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  constexpr auto ConstantOffset = offsetof(LocalIntegrationData, nApNm1);
  constexpr auto DissipationOffset = offsetof(NeighboringIntegrationData, nAmNm1);
  static_assert(ConstantOffset % sizeof(real) == 0 && DissipationOffset % sizeof(real) == 0,
                "A face's pair of matrices is not aligned to the real size.");

  constexpr std::array<FaceKinds, 2> BoundaryFaceKeys{FaceKinds::FreeSurface, FaceKinds::Outflow};

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany(*FaceRelations::Count, [&](void* stream, std::size_t faceRelation) {
      const ConditionalKey key(*KernelNames::NeighborFlux, *FaceKinds::Regular, face, faceRelation);
      if (table.find(key) == table.end()) {
        return;
      }
      auto& entry = table[key];

      kernel::gpu_damageLocalFlux local = deviceLocalFlux_;
      kernel::gpu_damageNeighborFlux neighbor = deviceNeighborFlux_;

      const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
      auto** dofs = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      const auto** own =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Integrals))->getDeviceDataPtr());
      const auto** other =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      const auto** localData = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      const auto** neighborData = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::NeighborIntegrationData))->getDeviceDataPtr());

      const auto bind = [&](auto& krnl) {
        krnl.numElements = numElements;
        krnl.Q = dofs;
        krnl.I = own;
        krnl.INeighbor = other;
        krnl.fluxConstant = localData;
        krnl.extraOffset_fluxConstant =
            (ConstantOffset / sizeof(real)) + face * tensor::fluxConstant::size();
        krnl.fluxDissipation = neighborData;
        krnl.extraOffset_fluxDissipation =
            (DissipationOffset / sizeof(real)) + face * tensor::fluxDissipation::size();
        krnl.streamPtr = stream;
      };

      const auto maxTmpMem = yateto::getMaxTmpMemRequired(local, neighbor);
      auto* tmpMem =
          reinterpret_cast<real*>(device_.api->allocMemAsync(maxTmpMem * numElements, stream));

      bind(local);
      local.linearAllocator.initialize(tmpMem);
      local.execute(face);

      bind(neighbor);
      neighbor.linearAllocator.initialize(tmpMem);
      (neighbor.*kernel::gpu_damageNeighborFlux::ExecutePtrs[faceRelation])();

      device_.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
    });
  }

  // A face without a neighbour, which is its own far state and its own wave
  // speed. Only the cell's own half of the flux runs for it -- there is no
  // second half to come -- and the ghost rule is already in the pair of
  // matrices, so the body is the local half of a regular face and nothing
  // else. It is keyed by the kind of the face rather than by a face
  // relation, because it has none.
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    runtime.envMany(BoundaryFaceKeys.size(), [&](void* stream, std::size_t kind) {
      const ConditionalKey key(
          *KernelNames::NeighborFlux, *BoundaryFaceKeys[kind], face, *FaceRelations::Any);
      if (table.find(key) == table.end()) {
        return;
      }
      auto& entry = table[key];

      kernel::gpu_damageLocalFlux local = deviceLocalFlux_;

      const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
      local.numElements = numElements;
      local.Q = (entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr();
      local.I =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Integrals))->getDeviceDataPtr());
      local.INeighbor =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr());
      local.fluxConstant = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      local.extraOffset_fluxConstant =
          (ConstantOffset / sizeof(real)) + face * tensor::fluxConstant::size();
      local.fluxDissipation = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::NeighborIntegrationData))->getDeviceDataPtr());
      local.extraOffset_fluxDissipation =
          (DissipationOffset / sizeof(real)) + face * tensor::fluxDissipation::size();
      local.streamPtr = stream;

      auto* tmpMem = reinterpret_cast<real*>(device_.api->allocMemAsync(
          kernel::gpu_damageLocalFlux::TmpMaxMemRequiredInBytes * numElements, stream));
      local.linearAllocator.initialize(tmpMem);
      local.execute(face);
      device_.api->freeMemAsync(reinterpret_cast<void*>(tmpMem), stream);
    });
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

std::pair<PerformanceEstimate, PerformanceEstimate>
    Neighbor::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
                      const std::array<CellDRMapping, Cell::NumFaces>& /*cellDrMapping*/) const {
  PerformanceEstimate neighbor;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (faceTypes[face] != FaceType::Regular) {
      continue;
    }
    // both traces of the face, the flux between them, and the lift back
    // both halves of the face's flux, each assembling the pair it applies
    neighbor += PerformanceEstimate::fromKernel<kernel::damageLocalFlux>(face);
    neighbor += PerformanceEstimate::fromKernel<kernel::damageNeighborFlux>(
        neighboringIndices[face][1], neighboringIndices[face][0], face);
  }

  // The rupture faces of a cell are the dynamic rupture kernel's, and it
  // reports them itself; what is counted here is the regular flux, which a
  // rupture face contributes nothing to.
  return {neighbor, PerformanceEstimate{}};
}

} // namespace seissol::kernels::solver::nonlinearck
