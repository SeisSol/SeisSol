// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_SURFACESTORAGE_H_
#define SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_SURFACESTORAGE_H_

#include "Alignment.h"
#include "Config.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Solver/BoundaryStorage.h"
namespace seissol {

struct SurfaceLTS {
  // held in the reals and the layout of the configuration of the layer
  template <typename Cfg>
  using FaceDisplacementArray = Real<Cfg>[tensor::faceDisplacement<Cfg>::size()];

  struct Side : public seissol::initializer::Variable<std::uint8_t> {};
  struct MeshId : public seissol::initializer::Variable<std::size_t> {};
  // the simulation copy of the cell of the face
  struct CopyId : public seissol::initializer::Variable<std::size_t> {};
  struct LocationFlag : public seissol::initializer::Variable<std::uint8_t> {};

  struct DisplacementDofs : public seissol::initializer::VariantVariable<FaceDisplacementArray> {};

  using Storage = initializer::Storage;
  using Layer = initializer::Layer;
  template <typename Cfg>
  using Ref = initializer::CellRef<Cfg>;
  using Backmap = initializer::StorageBackmap<1>;

  static void addTo(Storage& storage) {
    const seissol::initializer::LayerMask ghostMask(Ghost);
    storage.add<Side>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);
    storage.add<MeshId>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);
    storage.add<CopyId>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);
    storage.add<LocationFlag>(ghostMask, Alignment, initializer::AllocationMode::HostOnly);

    storage.add<DisplacementDofs>(ghostMask, PagesizeHeap, allocationModeBoundary());
  }
};

} // namespace seissol

#endif // SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_SURFACESTORAGE_H_
