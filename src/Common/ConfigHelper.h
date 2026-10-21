// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_CONFIGHELPER_H_
#define SEISSOL_SRC_COMMON_CONFIGHELPER_H_

#include "Common/Typedefs.h"

#include <string>
namespace seissol {

/// the name a basis has as a value of BASIS
std::string basisName(BasisType basis);

extern const std::string ConfigString;
extern const std::string ConfigDescriptor;

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_CONFIGHELPER_H_
