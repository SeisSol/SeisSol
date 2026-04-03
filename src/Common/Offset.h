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

// NOLINTBEGIN (bugprone-macro-parentheses)

// IMPORTANT: these macros use member names directly; and thus cannot be safeguarded by parentheses.

#define SEISSOL_STRINGIFY_INTERNAL(x) #x
#define SEISSOL_STRINGIFY(x) SEISSOL_STRINGIFY_INTERNAL(x)

#define SEISSOL_OFFSET(type, member) (offsetof(type, member) / sizeof(real))

#define SEISSOL_ARRAY_OFFSET(type, member, arrayidx)                                               \
  (SEISSOL_OFFSET(type, member[0]) +                                                               \
   (arrayidx) * (SEISSOL_OFFSET(type, member[1]) - SEISSOL_OFFSET(type, member[0])))

#define SEISSOL_OFFSET_ASSERT(type, member)                                                        \
  static_assert(offsetof(type, member) % sizeof(real) == 0,                                        \
                "Offset not compatible with the given real type. Type: " SEISSOL_STRINGIFY(        \
                    type) " . Member: " SEISSOL_STRINGIFY(member));

#define SEISSOL_ARRAY_OFFSET_ASSERT(type, member)                                                  \
  SEISSOL_OFFSET_ASSERT(type, member);                                                             \
  static_assert((offsetof(type, member[1]) - offsetof(type, member[0])) % sizeof(real) == 0,       \
                "Array offset not compatible with the given real type. Type: " SEISSOL_STRINGIFY(  \
                    type) " . Array member: " SEISSOL_STRINGIFY(member));

// NOLINTEND

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_OFFSET_H_
