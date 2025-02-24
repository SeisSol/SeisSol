// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_LOG2_H_
#define SEISSOL_SRC_INITIALIZER_TREE_LOG2_H_

namespace seissol::initializer {
template <unsigned N>
struct Log2 {
  static const unsigned Result = 1 + Log2<(N >> 1)>::Result;
};

template <>
struct Log2<1> {
  static const unsigned Result = 0;
};
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TREE_LOG2_H_
