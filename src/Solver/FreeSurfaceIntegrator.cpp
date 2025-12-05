// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "FreeSurfaceIntegrator.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <utils/logger.h>
#include <vector>

namespace seissol::solver {

FreeSurfaceIntegrator::FreeSurfaceIntegrator() = default;

FreeSurfaceIntegrator::~FreeSurfaceIntegrator() = default;

void FreeSurfaceIntegrator::initialize(unsigned /*maxRefinementDepth*/,
                                       LTS::Storage& ltsStorage,
                                       SurfaceLTS::Storage& surfaceStorage) {
  this->surfaceStorage = &surfaceStorage;

  m_enabled = true;

  logInfo() << "Initializing free surface integrator.";
  initializeSurfaceStorage(ltsStorage);
  logInfo() << "Initializing free surface integrator. Done.";
}

FreeSurfaceIntegrator::LocationFlag FreeSurfaceIntegrator::getLocationFlag(
    CellMaterialData materialData, FaceType faceType, unsigned int face) {
  if (initializer::isAcousticSideOfElasticAcousticInterface(materialData, face)) {
    return LocationFlag::Acoustic;
  } else if (initializer::isElasticSideOfElasticAcousticInterface(materialData, face)) {
    return LocationFlag::Elastic;
  } else if (faceType == FaceType::FreeSurface) {
    return LocationFlag::FreeSurface;
  } else if (faceType == FaceType::FreeSurfaceGravity) {
    return LocationFlag::FreeSurfaceWithGravity;
  } else {
    logError() << "Internal error in free surface integrator. Called for invalid cell.";
    std::abort(); // logError aborts, but this hides the compiler warning about missing return
                  // statement
  }
}

void FreeSurfaceIntegrator::initializeSurfaceStorage(LTS::Storage& ltsStorage) {
  const seissol::initializer::LayerMask ghostMask(Ghost);

  surfaceStorage->setLayerCount(ltsStorage.getColorMap());
  surfaceStorage->fixate();

  totalNumberOfFreeSurfaces = 0;
  for (auto [layer, surfaceLayer] :
       seissol::common::zip(ltsStorage.leaves(ghostMask), surfaceStorage->leaves(ghostMask))) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* cellMaterialData = layer.var<LTS::Material>();

    std::size_t numberOfFreeSurfaces = 0;
    std::size_t numberOfOutputFreeSurfaces = 0;
    const auto layerSize = layer.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)                                                          \
    reduction(+ : numberOfFreeSurfaces, numberOfOutputFreeSurfaces)
#endif // _OPENMP
    for (std::size_t cell = 0; cell < layerSize; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (initializer::requiresDisplacement(
                cellInformation[cell], cellMaterialData[cell], face)) {
          ++numberOfFreeSurfaces;

          if (secondaryInformation[cell].duplicate == 0) {
            ++numberOfOutputFreeSurfaces;
          }
        }
      }
    }
    surfaceLayer.setNumberOfCells(numberOfFreeSurfaces);
    totalNumberOfFreeSurfaces += numberOfOutputFreeSurfaces;
  }
  backmap.resize(totalNumberOfFreeSurfaces);

  surfaceStorage->allocateVariables();
  surfaceStorage->touchVariables();

  // NOTE: we store also for space storage duplicates here
  // thus, we need a non-duplicate lookup table (backmap)

  std::size_t surfaceCellOffset = 0; // Counts all surface cells of all layers
  std::size_t surfaceCellGlobal = 0;
  for (auto [layer, surfaceLayer] :
       seissol::common::zip(ltsStorage.leaves(ghostMask), surfaceStorage->leaves(ghostMask))) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    real(*dofs)[tensor::Q::size()] = layer.var<LTS::Dofs>();
    real*(*faceDisplacements)[4] = layer.var<LTS::FaceDisplacements>();
    real*(*faceDisplacementsDevice)[4] = layer.var<LTS::FaceDisplacementsDevice>();
    real** surfaceDofs = surfaceLayer.var<SurfaceLTS::Dofs>();
    auto* displacementDofs = surfaceLayer.var<SurfaceLTS::DisplacementDofs>();
    auto* displacementDofsDevice =
        surfaceLayer.var<SurfaceLTS::DisplacementDofs>(initializer::AllocationPlace::Device);
    auto* cellMaterialData = layer.var<LTS::Material>();
    auto* surfaceBoundaryMapping = surfaceLayer.var<SurfaceLTS::BoundaryMapping>();
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* locationFlagLayer = surfaceLayer.var<SurfaceLTS::LocationFlag>();

    auto* side = surfaceLayer.var<SurfaceLTS::Side>();
    auto* meshId = surfaceLayer.var<SurfaceLTS::MeshId>();
    auto* outputPosition = surfaceLayer.var<SurfaceLTS::OutputPosition>();
    std::size_t surfaceCell = 0;
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (initializer::requiresDisplacement(
                cellInformation[cell], cellMaterialData[cell], face)) {
          surfaceDofs[surfaceCell] = dofs[cell];

          // NOTE: assign LTS::Storage data here
          faceDisplacements[cell][face] = displacementDofs[surfaceCell];
          faceDisplacementsDevice[cell][face] = displacementDofsDevice[surfaceCell];

          side[surfaceCell] = face;
          meshId[surfaceCell] = secondaryInformation[cell].meshId;
          surfaceBoundaryMapping[surfaceCell] = &boundaryMapping[cell][face];
          locationFlagLayer[surfaceCell] = static_cast<std::uint8_t>(
              getLocationFlag(cellMaterialData[cell], cellInformation[cell].faceTypes[face], face));

          const auto globalId = secondaryInformation[cell].globalId * 4 + face;

          if (secondaryInformation[cell].duplicate == 0) {
            outputPosition[surfaceCell] = surfaceCellOffset;
            backmap[surfaceCellOffset] = surfaceCellGlobal;
            ++surfaceCellOffset;
          } else {
            outputPosition[surfaceCell] = std::numeric_limits<std::size_t>::max();
          }
          ++surfaceCell;
          ++surfaceCellGlobal;
        }
      }
    }
  }
}

} // namespace seissol::solver
