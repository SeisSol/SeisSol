// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_COMMON_CONFIGHELPER_H_
#define SEISSOL_SRC_COMMON_CONFIGHELPER_H_

#include <Config.h>
#include <array>

namespace seissol {

extern const std::array<ConfigVariant, std::variant_size_v<ConfigVariant>> ConfigVariantList;

extern const std::array<std::string, std::variant_size_v<ConfigVariant>> ConfigString;

} // namespace seissol
#endif // SEISSOL_SRC_COMMON_CONFIGHELPER_H_
