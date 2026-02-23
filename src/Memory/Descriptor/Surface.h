// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_

#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
namespace seissol {

struct SurfaceLTS {
  using FaceDisplacementType = real[tensor::faceDisplacement::size()];

  struct Dofs : public seissol::initializer::Variable<real*> {};
  struct Side : public seissol::initializer::Variable<std::uint8_t> {};
  struct MeshId : public seissol::initializer::Variable<std::size_t> {};
  struct BoundaryMapping : public seissol::initializer::Variable<CellBoundaryMapping*> {};
  struct LocationFlag : public seissol::initializer::Variable<std::uint8_t> {};

  struct DisplacementDofs : public seissol::initializer::Variable<FaceDisplacementType> {};

  struct SurfaceVarmap
      : public initializer::
            SpecificVarmap<Dofs, Side, MeshId, BoundaryMapping, LocationFlag, DisplacementDofs> {};

  using Storage = initializer::Storage<SurfaceVarmap>;
  using Layer = initializer::Layer<SurfaceVarmap>;
  using Ref = initializer::Layer<SurfaceVarmap>::CellRef;
  using Backmap = initializer::StorageBackmap<1>;

  static void addTo(Storage& storage) {
    const seissol::initializer::LayerMask ghostMask(Ghost);
    storage.add<Dofs>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    storage.add<Side>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    storage.add<MeshId>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    storage.add<BoundaryMapping>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    storage.add<LocationFlag>(ghostMask, 1, initializer::AllocationMode::HostOnly);

    storage.add<DisplacementDofs>(ghostMask, PagesizeHeap, allocationModeBoundary());
  }

  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          Storage& storage) {
    manager.registerData<DisplacementDofs>("displacementDofs", storage);
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
