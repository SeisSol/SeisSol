// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_PRECISION_H_
#define SEISSOL_SRC_KERNELS_PRECISION_H_

#include <Common/Real.h>
#include <Config.h>

// (real should be lower-case)

namespace seissol {
// NOLINTNEXTLINE
using real = RealT<Config::Precision>;
} // namespace seissol

namespace seissol {

enum class Precision { F32, F64 };

template <Precision P>
struct PrecisionWrapper {};

template <>
struct PrecisionWrapper<Precision::F32> {
  using Type = float;
};

template <>
struct PrecisionWrapper<Precision::F64> {
  using Type = double;
};

template <Precision P>
using PrecisionT = typename PrecisionWrapper<P>::Type;

} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_PRECISION_H_
