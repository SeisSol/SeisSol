// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_OFFSET_H_
#define SEISSOL_SRC_COMMON_OFFSET_H_

#include <cstddef>

namespace seissol {

#define SEISSOL_OFFSET(type, member) (offsetof(type, member) / sizeof(real))
#define SEISSOL_ARRAY_OFFSET(type, member, arrayidx)                                               \
  (SEISSOL_OFFSET(type, member[0]) +                                                               \
   (arrayidx) * (SEISSOL_OFFSET(type, member[1]) - SEISSOL_OFFSET(type, member[0])))

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_CONSTANTS_H_
