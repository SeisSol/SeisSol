// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_BOUNDARYHELPER_H_
#define SEISSOL_SRC_INITIALIZER_BOUNDARYHELPER_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

#include <cstddef>
#include <limits>

namespace seissol {

// TODO: make all constexpr with C++20.

inline bool isAcousticSideOfElasticAcousticInterface(CellMaterialData& material, std::size_t face) {
  constexpr auto Eps = std::numeric_limits<real>::epsilon();
  return material.neighbor[face]->getMuBar() > Eps && material.local->getMuBar() < Eps;
}
inline bool isElasticSideOfElasticAcousticInterface(CellMaterialData& material, std::size_t face) {
  constexpr auto Eps = std::numeric_limits<real>::epsilon();
  return material.local->getMuBar() > Eps && material.neighbor[face]->getMuBar() < Eps;
}

constexpr bool isAtElasticAcousticInterface(CellMaterialData& material, std::size_t face) {
  // We define the interface cells as all cells that are in the elastic domain but have a
  // neighbor with acoustic material.
  return material.local != nullptr && material.neighbor[face] != nullptr &&
         (isAcousticSideOfElasticAcousticInterface(material, face) ||
          isElasticSideOfElasticAcousticInterface(material, face));
}

constexpr bool requiresDisplacement(CellLocalInformation cellLocalInformation,
                                    CellMaterialData& material,
                                    std::size_t face) {
  const auto faceType = cellLocalInformation.faceTypes[face];
  return faceType == FaceType::FreeSurface || faceType == FaceType::FreeSurfaceGravity ||
         isAtElasticAcousticInterface(material, face);
}

constexpr bool requiresNodalFlux(FaceType f) {
  return (f == FaceType::FreeSurfaceGravity || f == FaceType::Dirichlet ||
          f == FaceType::Analytical);
}

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_BOUNDARYHELPER_H_
