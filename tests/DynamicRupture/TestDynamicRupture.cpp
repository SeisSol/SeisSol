// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "FrictionLaws/FrictionSolverCommon.t.h"
#include "Output/Geometry.t.h"
#include "Output/Variables.t.h"

namespace seissol::unit_test {
#define SEISSOL_CONFIGITER(cfg) TEST_CASE_TEMPLATE_INVOKE(configId2, cfg);
#include "ConfigInclude.h"

#define SEISSOL_CONFIGITER(cfg) TEST_CASE_TEMPLATE_INVOKE(configId3, cfg);
#include "ConfigInclude.h"

// TODO: update doctest first
// #define SEISSOL_CONFIGITER(cfg) TYPE_TO_STRING_AS(ConfigString[ConfigVariant(cfg()).index()],
// cfg); #include "ConfigInclude.h"

} // namespace seissol::unit_test
