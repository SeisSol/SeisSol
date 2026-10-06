// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FaultGeometry.h"

#include "DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"
#include "Geometry/CellGeometry.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshReader.h"
#include "Solver/MultipleSimulations.h"

#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::dr {

auto quadraturePoint(std::size_t point) -> geometry::FaceTransform::FaceVectorT {
  const auto points = init::quadpoints::view::create(init::quadpoints::Values);
  return {multisim::multisimTranspose(points, point, 0),
          multisim::multisimTranspose(points, point, 1)};
}

auto quadraturePointOf(std::size_t paddedPoint) -> std::optional<std::size_t> {
  const auto point = paddedPoint / multisim::NumSimulations;
  if (point >= misc::NumBoundaryGaussPoints) {
    return std::nullopt;
  }
  return point;
}

auto faultFramesAtPoints(std::size_t faultId, const geometry::MeshReader& mesh)
    -> std::vector<geometry::FaultFrame> {
  const auto& fault = mesh.getFault().at(faultId);
  const auto face = geometry::faultFaceTransformOf(faultId, mesh);
  std::vector<geometry::FaultFrame> frames;
  frames.reserve(misc::NumBoundaryGaussPoints);
  for (std::size_t point = 0; point < misc::NumBoundaryGaussPoints; ++point) {
    frames.push_back(geometry::faultFrameAt(*face, fault, quadraturePoint(point)));
  }
  return frames;
}

} // namespace seissol::dr
