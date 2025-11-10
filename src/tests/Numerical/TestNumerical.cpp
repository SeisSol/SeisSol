// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "BasisFunctions.t.h"
#include "Eigenvalues.t.h"
#include "Functions.t.h"
#include "GaussianNucleation.t.h"
#include "ODEInt.t.h"
#include "Quadrature.t.h"
#include "RegularizedYoffe.t.h"
#include "TimeBasis.t.h"
#include "Transformations.t.h"

namespace seissol::unit_test {
#define SEISSOL_CONFIGITER(cfg) TEST_CASE_TEMPLATE_INVOKE(configId5, cfg);
#include "ConfigInclude.h"

#define SEISSOL_CONFIGITER(cfg) TEST_CASE_TEMPLATE_INVOKE(configId7, cfg);
#include "ConfigInclude.h"

// TODO: update doctest first
// #define SEISSOL_CONFIGITER(cfg) TYPE_TO_STRING_AS(ConfigString[ConfigVariant(cfg()).index()],
// cfg); #include "ConfigInclude.h"

} // namespace seissol::unit_test
