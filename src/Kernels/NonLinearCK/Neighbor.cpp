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
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
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
  // read from it. Under global time stepping every cell keeps them, which is
  // the same condition local time stepping is turned off by.
  const real* own = data.get<LTS::StepIntegrals>();
  assert(own != nullptr);

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const auto faceType = info.faceTypes[face];
    const bool hasNeighbor = faceType == FaceType::Regular || faceType == FaceType::Periodic;

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
  logError() << "No GPU implementation provided";
}

std::pair<PerformanceEstimate, PerformanceEstimate>
    Neighbor::metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                      const std::array<std::array<uint8_t, 2>, Cell::NumFaces>& neighboringIndices,
                      const std::array<CellDRMapping, Cell::NumFaces>& /*cellDrMapping*/) const {
  PerformanceEstimate neighbor;

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (faceTypes[face] != FaceType::Regular && faceTypes[face] != FaceType::Periodic) {
      continue;
    }
    // both traces of the face, the flux between them, and the lift back
    // both halves of the face's flux, each assembling the pair it applies
    neighbor += PerformanceEstimate::fromKernel<kernel::damageLocalFlux>(face);
    neighbor += PerformanceEstimate::fromKernel<kernel::damageNeighborFlux>(
        neighboringIndices[face][1], neighboringIndices[face][0], face);
  }

  // dynamic rupture is not supported by a material whose traction is derived
  return {neighbor, PerformanceEstimate{}};
}

} // namespace seissol::kernels::solver::nonlinearck
