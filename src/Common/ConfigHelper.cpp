// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigHelper.h"

#include "Common/ConfigValue.h"
#include "Config.h"

#include <string>

namespace seissol {

const std::string ConfigString = configName(Config::Value);
const std::string ConfigDescriptor = describeConfig(Config::Value);

} // namespace seissol
