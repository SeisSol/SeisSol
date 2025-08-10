// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Attenuation.t.h"
#include "GodunovState.t.h"
#include <Common/ConfigHelper.h>

namespace seissol::unit_test {
#define SEISSOL_CONFIGITER(cfg) TEST_CASE_TEMPLATE_INVOKE(configId4, cfg);
#include "ConfigInclude.h"

// TODO: update doctest first
// #define SEISSOL_CONFIGITER(cfg) TYPE_TO_STRING_AS(ConfigString[ConfigVariant(cfg()).index()],
// cfg); #include "ConfigInclude.h"

} // namespace seissol::unit_test
