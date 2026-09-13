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
  logError() << "The face coupling of the nonlinear solver is not wired yet. Which shape it"
             << "takes follows from where the dissipation coefficient comes from: a coefficient"
             << "that the step reduces has to be read from both transported tensors here, and a"
             << "coefficient that the material bounds can be folded into a constant flux solver"
             << "per face at setup time instead.";
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
    neighbor += PerformanceEstimate::fromKernel<kernel::projectToFace>(face);
    neighbor += PerformanceEstimate::fromKernel<kernel::projectNeighborToFace>(
        neighboringIndices[face][1], neighboringIndices[face][0]);
    neighbor += PerformanceEstimate::fromKernel<kernel::damageRusanov>();
    neighbor += PerformanceEstimate::fromKernel<kernel::faceIntegral>(face);
  }

  // dynamic rupture is not supported by a material whose traction is derived
  return {neighbor, PerformanceEstimate{}};
}

} // namespace seissol::kernels::solver::nonlinearck
