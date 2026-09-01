// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Boundary.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Initializer/BoundaryHelper.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/Surface.h"
#include "Memory/Tree/Layer.h"
#include "Solver/FreeSurfaceIntegrator.h"

#include <cstddef>

namespace seissol::initializer::internal {

void initBoundaryStorage(Boundary::Storage& boundaryStorage, LTS::Storage& storage) {
  const LayerMask ghostMask(Ghost);

  boundaryStorage.setName("boundary");

  // Boundary face storage
  Boundary::addTo(boundaryStorage);
  boundaryStorage.setLayerCount(storage.getColorMap());
  boundaryStorage.fixate();

  // Iterate over layers of standard lts storage and face lts storage together.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(storage.leaves(ghostMask), boundaryStorage.leaves(ghostMask))) {
    const auto* cellInformation = layer.var<LTS::CellInformation>();

    std::size_t numberOfBoundaryFaces = 0;
    const auto layerSize = layer.size();

#pragma omp parallel for schedule(static) reduction(+ : numberOfBoundaryFaces)
    for (std::size_t cell = 0; cell < layerSize; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          ++numberOfBoundaryFaces;
        }
      }
    }
    boundaryLayer.setNumberOfCells(numberOfBoundaryFaces);
  }
  boundaryStorage.allocateVariables();
  boundaryStorage.touchVariables();

  // The boundary storage is now allocated, now we only need to map from cell lts
  // to face lts.
  // We do this by, once again, iterating over both storages at the same time.
  for (auto [layer, boundaryLayer] :
       seissol::common::zip(storage.leaves(ghostMask), boundaryStorage.leaves(ghostMask))) {
    const auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* boundaryMapping = layer.var<LTS::BoundaryMapping>();
    auto* boundaryMappingDevice = layer.var<LTS::BoundaryMappingDevice>();
    auto* faceInformation = boundaryLayer.var<Boundary::FaceInformation>(AllocationPlace::Host);
    auto* faceInformationDevice =
        boundaryLayer.var<Boundary::FaceInformation>(AllocationPlace::Device);

    std::size_t boundaryFace = 0;
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (requiresNodalFlux(cellInformation[cell].faceTypes[face])) {
          boundaryMapping[cell][face] = CellBoundaryMapping(faceInformation[boundaryFace]);
          boundaryMappingDevice[cell][face] =
              CellBoundaryMapping(faceInformationDevice[boundaryFace]);
          ++boundaryFace;
        } else {
          boundaryMapping[cell][face] = CellBoundaryMapping();
          boundaryMappingDevice[cell][face] = CellBoundaryMapping();
        }
      }
    }
  }
}

void initSurfaceStorage(SurfaceLTS::Storage& surfaceStorage,
                        LTS::Storage& storage,
                        solver::FreeSurfaceIntegrator& freeSurfaceIntegrator,
                        int refinement) {
  surfaceStorage.setName("surface");
  SurfaceLTS::addTo(surfaceStorage);

  // TODO: move freeSurfaceIntegrator initialization here, once separated from the IO (cf. #1180).
  freeSurfaceIntegrator.initialize(refinement, storage, surfaceStorage);
}

} // namespace seissol::initializer::internal
