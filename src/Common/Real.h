// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_REAL_H_
#define SEISSOL_SRC_COMMON_REAL_H_

#include <cstddef>

// cf. https://stackoverflow.com/a/70868019
#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__
#endif

#include <cfloat>

namespace seissol {

enum class RealType { F32, F64, F128 };

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

template <>
struct RealTypeWrapper<RealType::F128> {
  using Type = _Float128;
};

template <RealType P>
using RealT = typename RealTypeWrapper<P>::Type;

constexpr std::size_t sizeOfRealType(RealType type) {
  switch (type) {
  case seissol::RealType::F32:
    return 4;
  case seissol::RealType::F64:
    return 8;
  case seissol::RealType::F128:
    return 16;
  default:
    throw;
  }
}

} // namespace seissol

#include <iostream>

namespace std {

inline ostream& operator<<(ostream& stream, const _Float128& val) {
  stream << static_cast<double>(val);
  return stream;
}
} // namespace std

#endif // SEISSOL_SRC_COMMON_REAL_H_
