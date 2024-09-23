// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

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
