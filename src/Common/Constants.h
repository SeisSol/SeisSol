// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_CONSTANTS_H_
#define SEISSOL_SRC_COMMON_CONSTANTS_H_

#include "Config.h"
#include <cstddef>

namespace seissol {
// TODO: remove these, once properly templated
constexpr std::size_t ConvergenceOrder = Config::ConvergenceOrder;

constexpr std::size_t PagesizeHeap = 2097152;
constexpr std::size_t PagesizeStack = 4096;

constexpr auto zeroLengthArrayHandler(std::size_t x) -> std::size_t { return x == 0 ? 1 : x; }

template <typename T, std::size_t N>
using NZArray = T[zeroLengthArrayHandler(N)];

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONSTANTS_H_
