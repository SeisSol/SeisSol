// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_

#include "Alignment.h"
#include "Common/Real.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace seissol {

struct SurfaceLTS {
  // held in the reals and the layout of the configuration of the layer
  template <typename Cfg>
  using FaceDisplacementArray = Real<Cfg>[tensor::faceDisplacement<Cfg>::size()];

  struct Side : public seissol::initializer::Variable<std::uint8_t> {};
  struct MeshId : public seissol::initializer::Variable<std::size_t> {};
  struct LocationFlag : public seissol::initializer::Variable<std::uint8_t> {};

  struct DisplacementDofs : public seissol::initializer::VariantVariable<FaceDisplacementArray> {};

  /// The state of the derived outputs of the free surface, per face (cf.
  /// initializer::DerivedStateLayout).
  struct DerivedState : public seissol::initializer::Variable<double> {};

  struct SurfaceVarmap
      : public initializer::
            SpecificVarmap<Side, MeshId, LocationFlag, DisplacementDofs, DerivedState> {};

  using Storage = initializer::Storage<SurfaceVarmap>;
  using Layer = initializer::Layer<SurfaceVarmap>;
  template <typename Cfg>
  using Ref = initializer::Layer<SurfaceVarmap>::CellRef<Cfg>;
  using Backmap = initializer::StorageBackmap<1>;

  /// `derivedState`: values per face of the state of the derived outputs.
  static void addTo(Storage& storage, std::size_t derivedState = 0) {
    const seissol::initializer::LayerMask ghostMask(Ghost);
    storage.add<Side>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);
    storage.add<MeshId>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);
    storage.add<LocationFlag>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);

    storage.add<DisplacementDofs>(ghostMask, PagesizeHeap, allocationModeBoundary());
    // on the host only, as LTS::DerivedState
    storage.add<DerivedState>(derivedState > 0 ? ghostMask
                                               : ghostMask | initializer::LayerMask(Copy) |
                                                     initializer::LayerMask(Interior),
                              Alignment,
                              initializer::AllocationMode::HostOnly,
                              false,
                              std::max<std::size_t>(derivedState, 1));
  }

  static void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                          Storage& storage) {
    manager.registerData<DisplacementDofs>("displacementDofs", storage);
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_SURFACE_H_
