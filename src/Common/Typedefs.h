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

enum class ViscoMode { None, QuantityExtension, AnelasticTensor };

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_TYPEDEFS_H_
