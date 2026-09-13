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
#include <cstddef>
#include <cstdint>
#include <utility>
#include <utils/logger.h>

namespace seissol::kernels::solver::nonlinearck {

void Neighbor::setGlobalData(const CompoundGlobalData& global) {
  projectToFace_.bindGlobals(*global.onHost);
  projectNeighborToFace_.bindGlobals(*global.onHost);
  rusanov_.bindGlobals(*global.onHost);
  faceIntegral_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceProjectToFace_.bindGlobals(*global.onDevice);
  deviceProjectNeighborToFace_.bindGlobals(*global.onDevice);
  deviceRusanov_.bindGlobals(*global.onDevice);
  deviceFaceIntegral_.bindGlobals(*global.onDevice);
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

  // The wave speed rides in the last column of a transported tensor, in its
  // constant mode. Reading it costs no kernel -- but it goes through the
  // generated view, because the leading dimension is padded for alignment and
  // the column does not sit where the shape alone would put it.
  constexpr auto DissipationColumn = tensor::I::Shape[1] - 1;
  const real lambdaLocal = init::I::view::create(own)(0, DissipationColumn);

  alignas(Alignment) real plusData[tensor::AplusT::size()];
  alignas(Alignment) real minusData[tensor::AminusT::size()];

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (info.faceTypes[face] != FaceType::Regular && info.faceTypes[face] != FaceType::Periodic) {
      logError() << "The nonlinear solver has no boundary conditions yet; face type"
                 << static_cast<int>(info.faceTypes[face]) << "cannot be handled.";
    }

    // Both halves of the flux differ only in the sign of their dissipation,
    // and what they share is the larger of the two wave speeds. Neither side
    // reads the other's material for it.
    kernel::damageFluxDissipation dissipation = dissipation_;
    dissipation.fluxConstant = data.get<LTS::LocalIntegration>().nApNm1[face];
    const real lambdaNeighbor = init::I::view::create(timeIntegrated[face])(0, DissipationColumn);
    dissipation.lambdaMax = std::max(lambdaLocal, lambdaNeighbor);
    dissipation.AplusT = plusData;
    dissipation.AminusT = minusData;
    dissipation.execute();

    kernel::damageLocalFlux local = localFlux_;
    local.Q = data.get<LTS::Dofs>();
    local.I = own;
    local.AplusT = plusData;
    local.execute(face);

    kernel::damageNeighborFlux neighbor = neighborFlux_;
    neighbor.Q = data.get<LTS::Dofs>();
    neighbor.I = timeIntegrated[face];
    neighbor.AminusT = minusData;
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
    // the dissipation correction of the face, then both halves of its flux
    neighbor += PerformanceEstimate::fromKernel<kernel::damageFluxDissipation>();
    neighbor += PerformanceEstimate::fromKernel<kernel::damageLocalFlux>(face);
    neighbor += PerformanceEstimate::fromKernel<kernel::damageNeighborFlux>(
        neighboringIndices[face][1], neighboringIndices[face][0], face);
  }

  // dynamic rupture is not supported by a material whose traction is derived
  return {neighbor, PerformanceEstimate{}};
}

} // namespace seissol::kernels::solver::nonlinearck
