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

  seissol::initializer::Variable<real*> dofs;
  seissol::initializer::Variable<std::uint8_t> side;
  seissol::initializer::Variable<std::size_t> meshId;
  seissol::initializer::Variable<std::size_t> outputPosition;
  seissol::initializer::Variable<CellBoundaryMapping*> boundaryMapping;
  seissol::initializer::Variable<std::uint8_t> locationFlag;

  seissol::initializer::Variable<FaceDisplacementType> displacementDofs;

  void addTo(seissol::initializer::LTSTree& surfaceLtsTree) {
    const seissol::initializer::LayerMask ghostMask(Ghost);
    surfaceLtsTree.add(dofs, ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add(side, ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add(meshId, ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add(outputPosition, ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add(boundaryMapping, ghostMask, 1, initializer::AllocationMode::HostOnly);
    surfaceLtsTree.add(locationFlag, ghostMask, 1, initializer::AllocationMode::HostOnly);

    surfaceLtsTree.add(
        displacementDofs, ghostMask, PagesizeHeap, initializer::allocationModeBoundary());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   initializer::LTSTree* tree) const {
    manager.registerData("displacementDofs", tree, displacementDofs);
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
