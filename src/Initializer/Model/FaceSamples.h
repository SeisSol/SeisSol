// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_FACESAMPLES_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_FACESAMPLES_H_

#include "Common/Constants.h"
#include "Equations/Datastructures.h" // IWYU pragma: keep
#include "Memory/Tree/Backmap.h"
#include "Model/OperatorLayout.h"

#include <cstddef>
#include <vector>

namespace seissol::initializer {

/**
 * The material where the faces of the cells need it: at the points the flux of a face is
 * integrated at, and at the points of the fault rule on either side of a fault face. Both are
 * queried directly, each point pulled slightly into the cell it is asked for (see InwardPull), so
 * that the two sides of a face see the material of their own cells there. Held only while the
 * operators of the faces are set up.
 */
struct FaceSamples {
  using MaterialT = seissol::model::MaterialT;

  /// The material at the points of the faces of every queried cell, cell by cell, face by face.
  std::vector<MaterialT> flux;

  /// Per layer of the storage and cell of the layer, ghost cells included: which cell of `flux`
  /// it is.
  std::vector<std::vector<std::size_t>> cellOf;

  /// Per fault face of the mesh, in the order of MeshReader::getFault(): the material at the
  /// points of the fault rule on its plus side, then on its minus side.
  std::vector<MaterialT> fault;

  /// How many points a side of a fault face is sampled at.
  std::size_t faultPoints{0};

  /// The material at the points of the face `side` of the cell at `position`.
  [[nodiscard]] const MaterialT* atFace(const StoragePosition& position, std::size_t side) const {
    const auto cell = cellOf.at(position.color).at(position.cell);
    return &flux.at((cell * Cell::NumFaces + side) * FluxFaceNodes);
  }

  /// The material at the points of the fault face `face` on its plus (0) or minus (1) side.
  [[nodiscard]] const MaterialT* atFault(std::size_t face, std::size_t side) const {
    return &fault.at((2 * face + side) * faultPoints);
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_FACESAMPLES_H_
