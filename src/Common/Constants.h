// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONSTANTS_H_
#define SEISSOL_SRC_COMMON_CONSTANTS_H_

#include "Alignment.h"
#include "Config.h"
#include <cstddef>

namespace seissol {
struct Cell {
  static constexpr std::size_t NumFaces = 4;
  static constexpr std::size_t NumVertices = 4;
  static constexpr std::size_t Dim = 3;
};

constexpr auto zeroLengthArrayHandler(std::size_t x) -> std::size_t { return x == 0 ? 1 : x; }
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONSTANTS_H_
