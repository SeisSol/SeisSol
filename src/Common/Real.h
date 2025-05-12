// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_REAL_H_
#define SEISSOL_SRC_COMMON_REAL_H_

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

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_REAL_H_
