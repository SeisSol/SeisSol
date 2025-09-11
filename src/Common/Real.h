// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_REAL_H_
#define SEISSOL_SRC_COMMON_REAL_H_

#include <cstddef>
namespace seissol {

enum class RealType { F32, F64 };

template <RealType P>
struct RealTypeWrapper {};

template <>
struct RealTypeWrapper<RealType::F32> {
  using Type = float;
};

template <>
struct RealTypeWrapper<RealType::F64> {
  using Type = double;
};

template <RealType P>
using RealT = typename RealTypeWrapper<P>::Type;

constexpr std::size_t sizeOfRealType(RealType type) {
  switch (type) {
  case seissol::RealType::F32:
    return 4;
  case seissol::RealType::F64:
    return 8;
  default:
    throw;
  }
}

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_REAL_H_
