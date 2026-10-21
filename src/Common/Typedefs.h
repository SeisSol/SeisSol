// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_TYPEDEFS_H_
#define SEISSOL_SRC_COMMON_TYPEDEFS_H_

namespace seissol {

enum class DRQuadRuleType { Stroud, Dunavant, WitherdenVincent };

/// The basis of the polynomials of a cell and of a face (BASIS), with the matrices SeisSol/matrices
/// writes for it. A checkpoint stores the value, so a basis keeps its number.
enum class BasisType { Dubiner = 0, DubinerNormalized = 1, Bernstein = 2 };

/// The scheme that advances a cell in time. Which one a material may use is
/// decided at configure time, not here.
enum class SolverType { LinearCK, LinearCKAnelastic, STP };

enum class BuildType { Cpu, Gpu };

enum class DeviceBackend { None, Cuda, Hip, Sycl };

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_TYPEDEFS_H_
