// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Neighbor.h"

#include "Common/Constants.h"
#include "Common/Marker.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
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
  logError() << "The face coupling of the nonlinear solver is not wired yet: the constant half"
             << "of the flux solver has to be built per cell and face at setup, which is the"
             << "one thing still missing before the three kernel calls below can be made.";
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
