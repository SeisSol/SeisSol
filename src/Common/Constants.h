// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONSTANTS_H_
#define SEISSOL_SRC_COMMON_CONSTANTS_H_

#include "Alignment.h"

#include <array>
#include <cstddef>

namespace seissol {

struct Cell {
  static constexpr std::size_t NumFaces = 4;
  static constexpr std::size_t NumVertices = 4;
  static constexpr std::size_t Dim = 3;

  /// the barycenter of the reference cell, in reference coordinates
  static constexpr std::array<double, Dim> ReferenceBarycenter{
      1.0 / (Dim + 1), 1.0 / (Dim + 1), 1.0 / (Dim + 1)};
};

/// a face of a Cell; one dimension and one vertex less
struct Face {
  static constexpr std::size_t NumVertices = Cell::NumVertices - 1;
  static constexpr std::size_t Dim = Cell::Dim - 1;

  /// the number of distinct orientations a face can be seen under from its two adjacent cells
  static constexpr std::size_t NumOrientations = NumVertices;

  /// the barycenter of the reference face, in face coordinates
  static constexpr std::array<double, Dim> ReferenceBarycenter{1.0 / (Dim + 1), 1.0 / (Dim + 1)};
};

constexpr auto zeroGuard(std::size_t x) -> std::size_t { return x == 0 ? 1 : x; }
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONSTANTS_H_
