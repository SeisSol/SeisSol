// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigHelper.h"

#include <Config.h>
#include <array>

namespace seissol {

const std::array<ConfigVariant, std::variant_size_v<ConfigVariant>> ConfigVariantList {
#define _H_(cfg) cfg(),
#include "ConfigInclude.h"
};

} // namespace seissol
