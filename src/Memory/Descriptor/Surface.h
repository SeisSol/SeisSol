// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_

#include <Initializer/Typedefs.h>
#include <Memory/Descriptor/Boundary.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
namespace seissol {

struct SurfaceLTS {
  using FaceDisplacementType = real[tensor::faceDisplacement::size()];

  struct Dofs : public seissol::initializer::Variable<real*> {};
  struct Side : public seissol::initializer::Variable<std::uint8_t> {};
  struct MeshId : public seissol::initializer::Variable<std::size_t> {};
  struct OutputPosition : public seissol::initializer::Variable<std::size_t> {};
  struct BoundaryMapping : public seissol::initializer::Variable<CellBoundaryMapping*> {};
  struct LocationFlag : public seissol::initializer::Variable<std::uint8_t> {};

  struct DisplacementDofs : public seissol::initializer::Variable<FaceDisplacementType> {};

  struct SurfaceVarmap : public initializer::SpecificVarmap<Dofs,
                                                            Side,
                                                            MeshId,
                                                            OutputPosition,
                                                            BoundaryMapping,
                                                            LocationFlag,
                                                            DisplacementDofs> {};

  using Tree = initializer::LTSTree<SurfaceVarmap>;
  using Layer = initializer::Layer<SurfaceVarmap>;
  using Ref = initializer::Layer<SurfaceVarmap>::CellRef;
  using Backmap = initializer::StorageBackmap<1>;

  static void addTo(Tree& surfaceLtsTree) {
    const seissol::initializer::LayerMask ghostMask(Ghost);
    surfaceLtsTree.add<Dofs>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add<Side>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add<MeshId>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add<OutputPosition>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add<BoundaryMapping>(ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add<LocationFlag>(ghostMask, 1, initializer::AllocationMode::HostOnly);

    surfaceLtsTree.add<DisplacementDofs>(ghostMask, PagesizeHeap, allocationModeBoundary());
  }

  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          Tree* tree) {
    manager.registerData<DisplacementDofs>("displacementDofs", tree);
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
